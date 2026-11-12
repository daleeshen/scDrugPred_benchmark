# run_scad_benchmark.py
import argparse
import torch
import numpy as np
import pandas as pd
import sklearn.preprocessing as sk
import os
import itertools
import re
from itertools import cycle
from torch.utils.data.sampler import WeightedRandomSampler
import warnings
from collections import Counter

warnings.filterwarnings('ignore')

# 强制使用 CPU 保证资源测试公平
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

from SCADmodules import FX, MTLP, Discriminator


def extract_best_auc_params(log_text):
    pattern = re.compile(
        r"Drug = ([\w-]+), extractor = (\w+), source average validation AUC of iter(\d+) = ([\d.]+), "
        r"hdim(\d+)_zdim(\d+)_lr([\d.]+)_epoch(\d+)_lambda1([\d.]+)_"
        r"dropouts([\d.]+)_([\d.]+)_([\d.]+)_mbS(\d+)_mbT(\d+)"
    )
    best_auc = float('-inf')
    best_params = {}
    for match in pattern.finditer(log_text):
        drug, extractor, iter_num, auc, hdim, zdim, lr, epoch, lambda1, d1, d2, d3, mbS, mbT = match.groups()
        auc = float(auc)
        if auc > best_auc:
            best_auc = auc
            best_params = {
                "hdim": int(hdim), "zdim": int(zdim), "epoch": int(epoch),
                "lambda1": float(lambda1), "mbS": int(mbS), "mbT": int(mbT)
            }
    return best_params


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--log_file', type=str, required=True)
    parser.add_argument('--x_src', type=str, required=True)
    parser.add_argument('--y_src', type=str, required=True)
    parser.add_argument('--x_tar_tr', type=str, required=True)
    parser.add_argument('--y_tar_tr', type=str, required=True)
    parser.add_argument('--x_tar_te', type=str, required=True)
    args = parser.parse_args()

    device = torch.device('cpu')

    # 1. 提取最优参数
    with open(args.log_file, 'r') as f:
        best_params = extract_best_auc_params(f.read())

    if not best_params:
        print("❌ 未能从 Log 中解析到参数！")
        return

    h_dim, z_dim, epoch, lam1, mbS, mbT = (
        best_params['hdim'], best_params['zdim'], best_params['epoch'],
        best_params['lambda1'], best_params['mbS'], best_params['mbT']
    )
    lr, dropout_rate = 1e-3, 0.5

    # 2. 读取数据 (无需 XValGDSC/YValGDSC)
    XTrainGDSC = pd.read_csv(args.x_src, sep='\t', index_col=0, decimal='.')
    YTrainGDSC = pd.read_csv(args.y_src, sep='\t', index_col=0, decimal='.')
    XTrainCells = pd.read_csv(args.x_tar_tr, sep='\t', index_col=0, decimal='.')
    YTrainCells = pd.read_csv(args.y_tar_tr, sep='\t', index_col=0, decimal='.')
    XTestCells = pd.read_csv(args.x_tar_te, sep='\t', index_col=0, decimal='.')

    IE_dim = XTrainGDSC.shape[1]

    # 3. 标准化对齐
    XTrainCombined = pd.concat([XTrainGDSC, XTrainCells])
    scalerTrain = sk.StandardScaler()
    scalerTrain.fit(XTrainCombined.values)

    XTrainGDSC_N = scalerTrain.transform(XTrainGDSC.values)
    XTrainCells_N = scalerTrain.transform(XTrainCells.values)
    XTestCells_N = scalerTrain.transform(XTestCells.values)

    # 4. 构建 DataLoader
    class_counts = Counter(YTrainGDSC.iloc[:, 0])
    weight = {0: 1. / class_counts[0], 1: 1. / class_counts[1]}
    samples_weight = np.array([weight[t] for t in YTrainGDSC.iloc[:, 0].values])
    samples_weight = torch.from_numpy(samples_weight).reshape(-1)
    sampler = WeightedRandomSampler(samples_weight.type('torch.DoubleTensor'), len(samples_weight), replacement=True)

    CDataset = torch.utils.data.TensorDataset(torch.FloatTensor(XTrainGDSC_N), torch.FloatTensor(YTrainGDSC.values))
    PDataset = torch.utils.data.TensorDataset(torch.FloatTensor(XTrainCells_N),
                                              torch.FloatTensor(YTrainCells.values.astype(int)))

    CLoader = torch.utils.data.DataLoader(dataset=CDataset, batch_size=mbS, shuffle=False, sampler=sampler)
    PLoader = torch.utils.data.DataLoader(dataset=PDataset, batch_size=mbT, shuffle=True)

    # 5. 初始化模型
    Gen = FX(dropout_rate, IE_dim, h_dim, z_dim).to(device)
    Map = MTLP(dropout_rate, h_dim, z_dim).to(device)
    DG = Discriminator(dropout_rate, h_dim, z_dim).to(device)

    optimizer = torch.optim.Adagrad(itertools.chain(Gen.parameters(), Map.parameters(), DG.parameters()), lr=lr)
    C_loss = torch.nn.BCELoss()

    print(f"🚀 Started SCAD Adaptation (Epochs: {epoch}, Cells: {len(XTrainCells)})")

    # 6. 对抗训练 (核心消耗)
    for it in range(epoch):
        Gen.train();
        Map.train();
        DG.train()
        for i, data in enumerate(zip(cycle(CLoader), PLoader)):
            xs, ys = data[0][0].to(device), data[0][1].view(-1, 1).to(device)
            xt = data[1][0].to(device)
            if xs.size()[0] < 5 or xt.size()[0] < 5: continue

            F_xs, F_xt = Gen(xs), Gen(xt)
            yhat_xs = Map(F_xs)
            closs = C_loss(yhat_xs, ys)

            Labels, Labelt = torch.ones(F_xs.size(0), 1), torch.zeros(F_xt.size(0), 1)
            Lst, Xst = torch.cat([Labels, Labelt], 0).to(device), torch.cat([F_xs, F_xt], 0).to(device)

            DG_loss = C_loss(DG(Xst), Lst)
            Loss = closs + lam1 * DG_loss

            optimizer.zero_grad()
            Loss.backward()
            optimizer.step()

    # 7. 测试集前向推理
    Gen.eval();
    Map.eval()
    with torch.no_grad():
        TXTestCells_N = torch.FloatTensor(XTestCells_N).to(device)
        F_xt_test = Gen(TXTestCells_N)
        yhatt_test = Map(F_xt_test)

    print("✅ Inference complete.")


if __name__ == "__main__":
    main()