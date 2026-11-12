import os.path
import os
import torch #
from dataloader import KfoldDataset, dataloader
import argparse
from tqdm import tqdm
from model import CellDSFormer, SVM
from transformers import BertTokenizer
from torch import optim
from tool import Accuracy_score, F1_score, AUROC_score, Recall_score, Precision_score, AP_score, AMI, ARI
import pickle
import random
import numpy as np
import logging
import pandas as pd
import datetime


# 设置日志记录配置
def set_project_logger(proj_name, log_folder):
    proj_log_name = os.path.join(log_folder, proj_name + '_log.txt')

    # 获取根日志记录器
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # 清除之前的处理器，避免重复添加
    if logger.hasHandlers():
        logger.handlers.clear()

    # 添加文件处理器
    file_handler = logging.FileHandler(proj_log_name)
    file_handler.setLevel(logging.INFO)

    # 添加控制台处理器
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)

    # 设置日志格式
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    stream_handler.setFormatter(formatter)

    # 将处理器添加到日志记录器
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

def prepare(proj_name, data_path):
    parser = argparse.ArgumentParser(description='AI4Bio')
    parser.add_argument('--ep_num', type=int, default=5, help='epoch number of training')
    parser.add_argument('--train_batch_size', type=int, default=24, help='')
    parser.add_argument('--test_batch_size', type=int, default=24, help='')
    parser.add_argument('--data_path', type=str, default='./cell_dt',
                        help='subdt is just for test, plz use cell_dt. when using subdt do not caculate the AUROC')  # cell_dt

    if torch.cuda.is_available():
        device = 'cuda:0'
        print('Using GPU')
    else:
        device = 'cpu'
        print('Using CPU')

    parser.add_argument('--device', type=str, default=device, help='gpu or cpu')

    parser.add_argument('--lr', type=int, default=0.0001, help='')
    parser.add_argument('--folds', type=int, default=5, help='')
    parser.add_argument('--seq_length', type=int, default=2048, help='')
    parser.add_argument('--d_model', type=int, default=256, help='')
    parser.add_argument('--num_heads', type=int, default=8, help='')
    parser.add_argument('--num_heads_gat', type=int, default=3, help='')
    parser.add_argument('--d_ff', type=int, default=1024, help='')
    parser.add_argument('--vocab_size', type=int, default=18628, help='')
    parser.add_argument('--init_node_f', type=int, default=101, help='')
    parser.add_argument('--node_f', type=int, default=256, help='')
    parser.add_argument('--name', type=str, default='model', help='specify model name')

    args = parser.parse_args()

    with open('Gdata.pickle', 'rb') as file:
        Gdata = pickle.load(file)

    source_path = data_path + '/' + proj_name
    trdt_list, tedt_list = KfoldDataset(source_path, args.folds)

    return args, trdt_list, tedt_list, Gdata,


def run(proj_name, data_path, save_path):
    args, trdt_list, tedt_list, Gdata = prepare(proj_name, data_path=data_path)

    Gdata = Gdata.to(args.device)
    for f in [4, 3, 1, 2, 0]:
        model = CellDSFormer(args)
        #model = SVM(args)

        loss_function = torch.nn.CrossEntropyLoss()
        model_parameters = [parameter for parameter in model.parameters() if parameter.requires_grad]
        optimizer = optim.AdamW(model_parameters, lr=args.lr)

        model = model.to(args.device)

        train_data_loader, test_data_loader = dataloader(current_fold=f, train_list=trdt_list, test_list=tedt_list,
                                                         tr_bs=8, te_bs=8)

        train_labels_all, test_labels_all = [], []
        for batch in train_data_loader:
            _, labels = batch
            # labels是one-hot形式，比如 [1,0] 或 [0,1]
            label_idx = torch.argmax(labels, dim=1).numpy()
            train_labels_all.extend(label_idx.tolist())

        for batch in test_data_loader:
            _, labels = batch
            label_idx = torch.argmax(labels, dim=1).numpy()
            test_labels_all.extend(label_idx.tolist())

        train_counts = pd.Series(train_labels_all).value_counts().to_dict()
        test_counts = pd.Series(test_labels_all).value_counts().to_dict()

        print(f"\n📊 Fold {f} label distribution:")
        print(f"Train labels -> {train_counts}")
        print(f"Test  labels -> {test_counts}\n")

        for epoch in range(args.ep_num):
            loss_sum = 0.0
            Acc = []
            F1 = []
            AUROC = []
            Precision = []
            Recall = []
            APscore = []
            Ami = []
            Ari = []

            with tqdm(train_data_loader, ncols=80, position=0, leave=True) as batches:
                for b in batches:  # sample
                    all_train_labels = []
                    input_ids, labels = b  # batch_size*seq_len

                    label_idx = torch.argmax(labels, dim=1).detach().cpu().numpy()
                    all_train_labels.extend(label_idx.tolist())

                    if torch.all(torch.logical_or(torch.all(labels == torch.tensor([1, 0])), torch.all(labels == torch.tensor([0, 1]))))==True:
                        continue
                    input_ids = input_ids.to(args.device)
                    labels = labels.to(args.device)
                    pred = model(input_ids, Gdata)
                    labels = labels.float()
                    pred = pred.float()
                    loss = loss_function(pred, labels)

                    optimizer.zero_grad()
                    loss.backward()
                    optimizer.step()

                    loss_sum += loss.item()
                    acc = Accuracy_score(pred, labels)
                    f1 = F1_score(pred, labels)
                    aur = AUROC_score(pred, labels)   # not for subdt
                    pre = Precision_score(pred, labels)
                    rcl = Recall_score(pred, labels)
                    aps = AP_score(pred, labels)
                    ami = AMI(pred, labels)
                    ari = ARI(pred, labels)


                    Acc.append(acc)
                    F1.append(f1)
                    AUROC.append(aur)
                    Precision.append(pre)
                    Recall.append(rcl)
                    APscore.append(aps)
                    Ami.append(ami)
                    Ari.append(ari)

                label_series = pd.Series(all_train_labels)
                label_counts = label_series.value_counts().to_dict()
                logging.info(f"📊 Training Epoch {epoch} Fold {f} Label distribution: {label_counts}")

                avg_acc = sum(Acc) / len(Acc) if len(Acc) > 0 else 0
                avg_auroc = sum(AUROC) / len(AUROC) if len(AUROC) > 0 else 0
                avg_precision = sum(Precision) / len(Precision) if len(AUROC) > 0 else 0
                recall = sum(Recall) / len(Recall) if len(Recall) > 0 else 0
                f1 = sum(F1) / len(F1) if len(F1) > 0 else 0
                ap_score = sum(APscore) / len(APscore) if len(APscore) > 0 else 0
                ami = sum(Ami) / len(Ami) if len(Ami) > 0 else 0
                ari = sum(Ari) / len(Ari) if len(Ari) > 0 else 0

                logging.info(
                    f'Training epoch: {epoch}, Current_fold: {f}, loss: {loss_sum}, Accuracy: {avg_acc}, '
                    f'AUROC: {avg_auroc}, Precision: {avg_precision}, '
                    f'Recall: {recall}, F1: {f1}, APscore: {ap_score}, '
                    f'Ami: {ami}, Ari: {ari}')

            loss_sum = 0.0
            Acc = []
            F1 = []
            AUROC = []
            Precision = []
            Recall = []
            with torch.no_grad():
                with tqdm(test_data_loader, ncols=80, position=0, leave=True) as batches:
                    all_test_labels = []

                    for b in batches:  # sample
                        input_ids, labels = b  # batch_size*seq_len

                        label_idx = torch.argmax(labels, dim=1).detach().cpu().numpy()
                        all_test_labels.extend(label_idx.tolist())

                        if torch.all(torch.logical_or(torch.all(labels == torch.tensor([1, 0])),
                                                      torch.all(labels == torch.tensor([0, 1])))) == True:
                            continue
                        input_ids = input_ids.to(args.device)
                        labels = labels.to(args.device)
                        pred = model(input_ids, Gdata)
                        labels = labels.float()
                        pred = pred.float()
                        loss = loss_function(pred, labels)

                        loss_sum += loss.item()
                        acc = Accuracy_score(pred, labels)
                        f1 = F1_score(pred, labels)
                        aur = AUROC_score(pred, labels)  # not for subdt
                        pre = Precision_score(pred, labels)
                        rcl = Recall_score(pred, labels)
                        aps = AP_score(pred, labels)
                        ami = AMI(pred, labels)
                        ari = ARI(pred, labels)


                        Acc.append(acc)
                        F1.append(f1)
                        AUROC.append(aur)
                        Precision.append(pre)
                        Recall.append(rcl)
                        APscore.append(aps)
                        Ami.append(ami)
                        Ari.append(ari)

                    label_series = pd.Series(all_test_labels)
                    label_counts = label_series.value_counts().to_dict()
                    logging.info(f"🧪 Testing Epoch {epoch} Fold {f} Label distribution: {label_counts}")

                    mean_ap = sum(APscore) / len(APscore)
                    baseline = np.mean(all_test_labels)  # 正样本比例
                    norm_ap = (mean_ap - baseline) / (1 - baseline) if baseline < 1 else 0
                    avg_acc = sum(Acc) / len(Acc) if len(Acc) > 0 else 0
                    avg_auroc = sum(AUROC) / len(AUROC) if len(AUROC) > 0 else 0
                    avg_precision = sum(Precision) / len(Precision) if len(AUROC) > 0 else 0
                    recall = sum(Recall) / len(Recall) if len(Recall) > 0 else 0
                    f1 = sum(F1) / len(F1) if len(F1) > 0 else 0
                    ap_score = sum(APscore) / len(APscore) if len(APscore) > 0 else 0
                    ami = sum(Ami) / len(Ami) if len(Ami) > 0 else 0
                    ari = sum(Ari) / len(Ari) if len(Ari) > 0 else 0
                    logging.info(f'Testing epoch: {epoch}, Current_fold: {f}, loss: {loss_sum}, '
                                 f'Accuracy: {avg_acc}, AUROC: {avg_auroc}, '
                                 f'Precision: {avg_precision}, Recall: {recall}, '
                                 f'F1: {f1}, APscore: {ap_score}, Normalized_AP: {norm_ap:.4f},'
                                 f'Ami: {ami}, Ari: {ari}')

        current_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        model_dir = save_path + '/model_save/'+proj_name + '/'
        if not os.path.exists(model_dir):
            os.makedirs(model_dir)
            print('model dir for ' + proj_name + ' has been created.')

        model_save_path = os.path.join(model_dir, f"{args.name}_{current_time}.ckpt")
        torch.save(model.state_dict(), model_save_path)


def set_seed(seed=42):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)  # 为所有CUDA设备设置随机种子
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

# 遍历数据集
root_path = '/media/shenli/Expansion/scDrug_benchmark/DrugFormer'
k_ratio_list = ['1_3', '1_10', '1_30', '1_100']
dataset_info = pd.read_csv('./dataset_drug_single.csv')

for round_idx in range(1, 4):
    for k_ratio in k_ratio_list:

        # ===== 跳过规则 =====
        # round 1 整体跳过
        if round_idx == 1:
            print(f'⏩ 跳过 round {round_idx}（全部 ratio）')
            continue

        # round 2 的 1_3、1_10、1_30 整体跳过
        if round_idx == 2 and k_ratio in {'1_3', '1_10', '1_30'}:
            print(f'⏩ 跳过 round {round_idx}, ratio {k_ratio}')
            continue
        # ===== 跳过规则结束 =====

        log_folder = os.path.join(root_path, f'{round_idx}_cell_dt_{k_ratio}_logs')
        save_path = os.path.join(root_path, f'{round_idx}_cell_dt_{k_ratio}_models')
        data_path = os.path.join(root_path, f'{round_idx}_{k_ratio}_cell_dt')

        # 创建日志文件夹和模型保存文件夹（如果不存在）
        if not os.path.exists(log_folder):
            os.makedirs(log_folder, exist_ok=True)
            print(f'Log folder created: {log_folder}')

        if not os.path.exists(save_path):
            os.makedirs(save_path, exist_ok=True)
            print(f'Model save folder created: {save_path}')

        # 数据集选择规则
        if round_idx == 2 and k_ratio == '1_100':
            # 从第10个数据集开始（索引9）
            run_list = dataset_info.iloc[9:, 0]
            print(f'⚙️ 从第9个数据集开始: round={round_idx}, ratio={k_ratio}')
        else:
            # round 3 正常；round 2 其它未被跳过的 ratio 也正常（如果有）
            run_list = dataset_info.iloc[:, 0]

        for ID in run_list:
            project_name = ID
            set_project_logger(project_name, log_folder=log_folder)
            print('开始处理数据集：', str(project_name))
            set_seed(seed=42)
            try:
                run(project_name, data_path=data_path, save_path=save_path)
                print(f'✅ 数据集 {project_name} 运行完成\n')
            except FileNotFoundError:
                print(f'⚠️ [FileNotFoundError] 数据集 {project_name} 缺失，跳过')
                logging.warning(f'[FileNotFoundError] Dataset folder not found: {data_path}')
                continue
            except Exception as e:
                # 捕获所有异常，打印具体错误但不中断循环
                print(f'❌ [Error] 数据集 {project_name} 运行失败：{type(e).__name__} - {e}')
                logging.error(f'[Error] Dataset {project_name} failed: {e}', exc_info=True)
                continue
print(dataset_info.iloc[:, 0])
'''
for data_path in data_path_list:
    if data_path == './cell_dt_2nd_round':
        log_folder = './2nd_logs/'
        save_path = './2nd_models'
    else:
        log_folder = './3rd_logs/'
        save_path = './3rd_models'
    for ID in dataset_info.iloc[:, 0]:
        proj_name = ID
        proj_log_name = log_folder + proj_name + '_log.txt'
        set_project_logger(proj_name)

        print('开始处理数据集：', str(proj_name))
        set_seed(seed=42)
        run(proj_name, data_path=data_path, save_path=save_path)


'''
