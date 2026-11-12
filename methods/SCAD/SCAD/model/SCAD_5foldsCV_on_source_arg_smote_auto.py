import argparse
import torch
import numpy as np
import pandas as pd
import os
import random
import itertools
from collections import Counter
from itertools import cycle
from sklearn.metrics import *
import sklearn.preprocessing as sk
from torch.optim.lr_scheduler import ReduceLROnPlateau
from torch.utils.data.sampler import WeightedRandomSampler
from imblearn.over_sampling import SMOTE
from imblearn.pipeline import Pipeline
from imblearn.under_sampling import RandomUnderSampler
from SCADmodules import *

def get_args():
    parser = argparse.ArgumentParser(description="SCAD Model Training with SMOTE Hyperparameter Search")
    parser.add_argument('--id', type=str, required=True, help='Dataset ID (e.g., GSE140440)')
    parser.add_argument('--drug', type=str, required=True, help='Drug name (e.g., Docetaxel)')
    parser.add_argument('--solid_tumor', type=int, choices=[0, 1], required=True, help='1 for solid tumor, 0 for non-solid tumor')
    parser.add_argument('--base_path', type=str, required=True, help='Base directory containing processed data and results')
    parser.add_argument('--gpu', type=int, default=1, choices=[0, 1], help='1 for CUDA, 0 for CPU')
    parser.add_argument('--max_iter', type=int, default=20, help='Number of random hyperparameter search iterations')
    return parser.parse_args()

def roc_auc_score_trainval(y_true, y_predicted):
    if len(np.unique(y_true)) == 1:
        return accuracy_score(y_true, np.rint(y_predicted))
    return roc_auc_score(y_true, y_predicted)

def reset_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True

def main():
    args = get_args()

    device = "cuda" if args.gpu == 1 and torch.cuda.is_available() else "cpu"
    print(f"Using device: {device}")

    folder_name = f"{args.id}_{args.drug}"
    Drug = args.drug
    extractor = 'FX'
    init_seed = 42

    ls_splits = ['split1', 'split2', 'split3', 'split4', 'split5']
    ls_mb_size = [{'mbS': 8, 'mbT': 8}, {'mbS': 32, 'mbT': 32}]
    ls_h_dim = [1024, 512]
    ls_z_dim = [256, 128]
    ls_lr = [0.001]
    ls_epoch = [10, 20, 30, 40, 60, 80, 100]
    ls_lam = [0.2, 0.6, 0.8, 1.0, 2, 3, 4, 5]
    ls_dropout_gen = [0.5]
    ls_dropout_mtl = ls_dropout_gen
    ls_dropout_dg = ls_dropout_gen

    genesets = ['_norm', '_ppi', '_tp4k']

    for gs in genesets:
        # The original script runs each geneset in a separate process.
        reset_seed(init_seed)

        print("================================================================================")
        print(f"SMOTE hyperparameter search: {folder_name} | geneset={gs}")
        print("================================================================================")

        results_dir = os.path.join(args.base_path, "results_smote" + gs)
        summary_path = os.path.join(
            results_dir,
            'train_eval_on_source_hypers_summary_5fds_smote_fixlr.txt'
        )

        mode = "fine_tune/" + extractor
        save_results_to = os.path.join(results_dir, folder_name, mode)
        save_trace_to = os.path.join(save_results_to, "trace")

        os.makedirs(os.path.join(save_results_to, 'model'), exist_ok=True)
        os.makedirs(os.path.join(save_results_to, 'test'), exist_ok=True)

        load_data_from = os.path.join(
            args.base_path,
            'data',
            'split_norm',
            folder_name,
            'stratified' + gs
        )

        if gs == '_norm':
            source_dir = 'source_5_folds_norm'
            target_dir = 'target_5_folds_norm'
            source_file_name = f'Source_exprs_resp_z.{Drug}.tsv'
        else:
            source_dir = 'source_5_folds_' + gs
            target_dir = 'target_5_folds_' + gs
            source_file_name = f'Source_exprs_resp_z.{Drug}{gs}.tsv'

        source_data_path = os.path.join(
            args.base_path,
            'data',
            'split_norm',
            folder_name,
            source_file_name
        )

        YGDSC = pd.read_csv(
            source_data_path,
            sep='\t',
            index_col=0,
            decimal='.'
        )

        class_sample_count = np.array([
            Counter(YGDSC['response'])[0] / len(YGDSC['response']),
            Counter(YGDSC['response'])[1] / len(YGDSC['response'])
        ])

        # Same SMOTE pipeline as the original implementation.
        over = SMOTE(sampling_strategy=0.5)
        under = RandomUnderSampler(sampling_strategy=0.5)
        pipeline = Pipeline(steps=[('o', over), ('u', under)])

        for mbsize in ls_mb_size:
            mbS = mbsize['mbS']
            mbT = mbsize['mbT']

            for iters in range(args.max_iter):
                print(f"\nITERATION # {iters} | mbS={mbS} | mbT={mbT}")

                # Preserve the original random-choice sequence.
                hdm = random.choice(ls_h_dim)
                zdm = random.choice(ls_z_dim)
                lrs = random.choice(ls_lr)
                epch = random.choice(ls_epoch)
                lambd1 = random.choice(ls_lam)

                drop_gen = random.choice(ls_dropout_gen)
                drop_mtl = random.choice(ls_dropout_mtl)
                drop_dg = random.choice(ls_dropout_dg)

                h_dim = hdm
                z_dim = zdm
                lr = lrs
                epoch = epch
                lam1 = lambd1

                dropout_gen = drop_gen
                dropout_mtl = drop_gen
                dropout_dg = drop_gen

                print(
                    f"h_dim={h_dim}, z_dim={z_dim}, lr={lr}, "
                    f"epoch={epoch}, lambda1={lam1}, "
                    f"mbS={mbS}, mbT={mbT}"
                )

                hyper_para_record = (
                    f"hdim{h_dim}_zdim{z_dim}_lr{lr}_epoch{epoch}_"
                    f"lambda1{lam1}_dropouts{dropout_gen}_"
                    f"{dropout_mtl}_{dropout_dg}_mbS{mbS}_mbT{mbT}"
                )

                batch_sizes = f"mbS{mbS}_mbT{mbT}"
                model_params = (
                    f"hdim{h_dim}_zdim{z_dim}_lr{lr}_epoch{epoch}_"
                    f"lamb1{lam1}_dropouts{dropout_gen}_"
                    f"{dropout_mtl}_{dropout_dg}_{batch_sizes}"
                )

                AUCval_source_splits_total = []

                for split in ls_splits:
                    AUCval_source_1split = []

                    print(f"Reading data for {split}...")

                    XTrainGDSC = pd.read_csv(
                        os.path.join(load_data_from, source_dir, split, 'X_train_source.tsv'),
                        sep='\t',
                        index_col=0,
                        decimal='.'
                    )
                    YTrainGDSC = pd.read_csv(
                        os.path.join(load_data_from, source_dir, split, 'Y_train_source.tsv'),
                        sep='\t',
                        index_col=0,
                        decimal='.'
                    )

                    # SMOTE sampling on source training data.
                    XTrainGDSC, YTrainGDSC = pipeline.fit_resample(
                        XTrainGDSC,
                        YTrainGDSC
                    )

                    XValGDSC = pd.read_csv(
                        os.path.join(load_data_from, source_dir, split, 'X_val_source.tsv'),
                        sep='\t',
                        index_col=0,
                        decimal='.'
                    )
                    YValGDSC = pd.read_csv(
                        os.path.join(load_data_from, source_dir, split, 'Y_val_source.tsv'),
                        sep='\t',
                        index_col=0,
                        decimal='.'
                    )
                    XTrainCells = pd.read_csv(
                        os.path.join(load_data_from, target_dir, split, 'X_train_target.tsv'),
                        sep='\t',
                        index_col=0,
                        decimal='.'
                    )
                    YTrainCells = pd.read_csv(
                        os.path.join(load_data_from, target_dir, split, 'Y_train_target.tsv'),
                        sep='\t',
                        index_col=0,
                        decimal='.'
                    )
                    XTestCells = pd.read_csv(
                        os.path.join(load_data_from, target_dir, split, 'X_test_target.tsv'),
                        sep='\t',
                        index_col=0,
                        decimal='.'
                    )
                    YTestCells = pd.read_csv(
                        os.path.join(load_data_from, target_dir, split, 'Y_test_target.tsv'),
                        sep='\t',
                        index_col=0,
                        decimal='.'
                    )

                    dir_name = os.path.join(
                        save_trace_to,
                        batch_sizes,
                        model_params
                    )
                    os.makedirs(dir_name, exist_ok=True)

                    XTrainCombined = pd.concat([
                        XTrainGDSC,
                        XTrainCells
                    ])

                    scalerTrain = sk.StandardScaler()
                    scalerTrain.fit(XTrainCombined.values)

                    XTrainGDSC_N = scalerTrain.transform(XTrainGDSC.values)
                    XValGDSC_N = scalerTrain.transform(XValGDSC.values)
                    XTrainCells_N = scalerTrain.transform(XTrainCells.values)
                    XTestCells_N = scalerTrain.transform(XTestCells.values)

                    TXValGDSC_N = torch.FloatTensor(
                        XValGDSC_N
                    ).to(device)
                    TYValGDSC = torch.FloatTensor(
                        YValGDSC.values
                    ).to(device)

                    # The original code constructs this sampler but does not
                    # actually pass it to the source DataLoader.
                    weight = 1. / class_sample_count
                    samples_weight = np.array([
                        weight[t] for t in YTrainGDSC.values
                    ])
                    samples_weight = torch.from_numpy(
                        samples_weight
                    ).reshape(-1)

                    sampler = WeightedRandomSampler(
                        samples_weight.type('torch.DoubleTensor'),
                        len(samples_weight),
                        replacement=True
                    )

                    CDataset = torch.utils.data.TensorDataset(
                        torch.FloatTensor(XTrainGDSC_N),
                        torch.FloatTensor(YTrainGDSC.values)
                    )

                    PDataset = torch.utils.data.TensorDataset(
                        torch.FloatTensor(XTrainCells_N),
                        torch.FloatTensor(YTrainCells.values.astype(int))
                    )

                    # Preserve original behavior: sampler is not used here.
                    CLoader = torch.utils.data.DataLoader(
                        dataset=CDataset,
                        batch_size=mbS,
                        shuffle=True
                    )

                    PLoader = torch.utils.data.DataLoader(
                        dataset=PDataset,
                        batch_size=mbT,
                        shuffle=True
                    )

                    _, IE_dim = XTrainGDSC_N.shape

                    if extractor == 'FX' or extractor == 'FX_MLP':
                        Gen = FX(
                            dropout_gen,
                            IE_dim,
                            h_dim,
                            z_dim
                        )
                    else:
                        raise NotImplementedError(
                            f'{extractor} is not a valid extractor'
                        )

                    Map = MTLP(
                        dropout_mtl,
                        h_dim,
                        z_dim
                    )
                    DG = Discriminator(
                        dropout_dg,
                        h_dim,
                        z_dim
                    )

                    Gen.to(device)
                    Map.to(device)
                    DG.to(device)

                    optimizer = torch.optim.Adagrad(
                        itertools.chain(
                            Gen.parameters(),
                            Map.parameters(),
                            DG.parameters()
                        ),
                        lr=lr
                    )

                    # Present in the original implementation but never stepped.
                    scheduler = ReduceLROnPlateau(
                        optimizer,
                        patience=5
                    )

                    C_loss = torch.nn.BCELoss()

                    for it in range(epoch):
                        AUCvals_source = []

                        for data in zip(
                            cycle(CLoader),
                            PLoader
                        ):
                            DataS = data[0]
                            DataT = data[1]

                            xs = DataS[0].to(device)
                            ys = DataS[1].view(-1, 1).to(device)
                            xt = DataT[0].to(device)
                            yt = DataT[1].view(-1, 1).to(device)

                            if xs.size()[0] < 5 or xt.size()[0] < 5:
                                continue

                            Gen.train()
                            Map.train()
                            DG.train()

                            if extractor == 'FX' or extractor == 'FX_MLP':
                                F_xs = Gen(xs)
                                F_xt = Gen(xt)
                            else:
                                raise NotImplementedError(
                                    f'{extractor} is not a valid extractor'
                                )

                            yhat_xs = Map(F_xs)
                            closs = C_loss(
                                yhat_xs,
                                ys
                            )

                            Labels = torch.ones(
                                F_xs.size(0),
                                1
                            )
                            Labelt = torch.zeros(
                                F_xt.size(0),
                                1
                            )

                            Lst = torch.cat(
                                [Labels, Labelt],
                                0
                            ).to(device)

                            Xst = torch.cat(
                                [F_xs, F_xt],
                                0
                            ).to(device)

                            yhat_DG = DG(Xst)
                            DG_loss = C_loss(
                                yhat_DG,
                                Lst
                            )

                            loss1 = closs
                            loss2 = lam1 * DG_loss
                            Loss = loss1 + loss2

                            optimizer.zero_grad()
                            Loss.backward()
                            optimizer.step()

                        with torch.no_grad():
                            Gen.eval()
                            Gen.to(device)
                            Map.eval()
                            Map.to(device)
                            DG.eval()
                            DG.to(device)

                            if extractor == 'FX' or extractor == 'FX_MLP':
                                F_xs_val = Gen(
                                    TXValGDSC_N
                                )
                            else:
                                raise NotImplementedError(
                                    f'{extractor} is not a valid extractor'
                                )

                            yhats_val = Map(
                                F_xs_val
                            )

                            ys_true_val = TYValGDSC.view(
                                -1,
                                1
                            ).cpu()

                            yhats_val = yhats_val.cpu()

                            AUC_source_val = roc_auc_score_trainval(
                                ys_true_val.detach().numpy(),
                                yhats_val.detach().numpy()
                            )

                            AUCvals_source.append(
                                AUC_source_val
                            )
                            AUCval_source_splits_total.append(
                                AUC_source_val
                            )
                            AUCval_source_1split.append(
                                AUC_source_val
                            )

                    avgAUCval_source_1split = np.mean(
                        np.array(AUCval_source_1split)
                    )

                    print(
                        f"{split}: source average validation AUC = "
                        f"{avgAUCval_source_1split:.6f}"
                    )

                avgAUCval_source = np.mean(
                    np.array(AUCval_source_splits_total)
                )

                with open(summary_path, 'a') as f:
                    f.write(
                        f"Drug={Drug}, Extractor={extractor}, "
                        f"Iter={iters}, Avg_Val_AUC={avgAUCval_source:.6f}, "
                        f"Params={hyper_para_record}\n"
                    )

                print(
                    f"Iteration {iters}: average source validation AUC = "
                    f"{avgAUCval_source:.6f}"
                )

if __name__ == "__main__":
    main()