import os
import torch
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


def set_project_logger(proj_name, log_folder):
    proj_log_name = os.path.join(log_folder, proj_name + '_log.txt')

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    if logger.hasHandlers():
        logger.handlers.clear()

    file_handler = logging.FileHandler(proj_log_name)
    file_handler.setLevel(logging.INFO)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    stream_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)


def prepare(proj_name, data_path):
    parser = argparse.ArgumentParser(description='AI4Bio')
    parser.add_argument('--ep_num', type=int, default=3, help='epoch number of training')
    parser.add_argument('--train_batch_size', type=int, default=12, help='')
    parser.add_argument('--test_batch_size', type=int, default=24, help='')
    parser.add_argument('--data_path', type=str, default='./subdt',
                        help='subdt is just for test, plz use cell_dt. when using subdt do not caculate the AUROC')
    parser.add_argument('--device', type=str, default='cuda:2' if torch.cuda.is_available() else 'cpu', help='')
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

    args, unknown = parser.parse_known_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    gdata_file = os.path.join(script_dir, 'Gdata.pickle')
    with open(gdata_file, 'rb') as file:
        Gdata = pickle.load(file)

    source_path = os.path.join(data_path, proj_name)
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
                                                         tr_bs=args.train_batch_size, te_bs=args.test_batch_size)

        for epoch in range(args.ep_num):
            loss_sum = 0
            Acc = []
            F1 = []
            AUROC = []
            Precision = []
            Recall = []
            APscore = []
            Ami = []
            Ari = []

            with tqdm(train_data_loader, ncols=80, position=0, leave=True) as batches:
                for b in batches:
                    input_ids, labels = b

                    if torch.all(torch.logical_or(torch.all(labels == torch.tensor([1, 0])),
                                                  torch.all(labels == torch.tensor([0, 1])))) == True:
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

                    loss_sum = loss_sum + loss
                    acc = Accuracy_score(pred, labels)
                    f1 = F1_score(pred, labels)
                    aur = AUROC_score(pred, labels)
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

                loss_value = loss_sum.item() if torch.is_tensor(loss_sum) else loss_sum
                logging.info(
                    f'Training epoch: {epoch}, Current_fold: {f}, loss: {loss_value}, Accuracy: {sum(Acc) / len(Acc)}, '
                    f'AUROC: {sum(AUROC) / len(AUROC)}, Precision: {sum(Precision) / len(Precision)}, '
                    f'Recall: {sum(Recall) / len(Recall)}, F1: {sum(F1) / len(F1)}, APscore: {sum(APscore) / len(APscore)}, '
                    f'Ami: {sum(Ami) / len(Ami)}, Ari: {sum(Ari) / len(Ari)}')

            loss_sum = 0
            Acc = []
            F1 = []
            AUROC = []
            Precision = []
            Recall = []
            all_test_labels = []

            with torch.no_grad():
                with tqdm(test_data_loader, ncols=80, position=0, leave=True) as batches:
                    for b in batches:
                        input_ids, labels = b

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

                        loss_sum = loss_sum + loss
                        acc = Accuracy_score(pred, labels)
                        f1 = F1_score(pred, labels)
                        aur = AUROC_score(pred, labels)
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

                    mean_ap = sum(APscore) / len(APscore)
                    baseline = np.mean(all_test_labels)
                    norm_ap = (mean_ap - baseline) / (1 - baseline) if baseline < 1 else 0
                    loss_value = loss_sum.item() if torch.is_tensor(loss_sum) else loss_sum

                    logging.info(
                        f'Testing epoch: {epoch}, Current_fold: {f}, loss: {loss_value}, '
                        f'Accuracy: {sum(Acc) / len(Acc)}, AUROC: {sum(AUROC) / len(AUROC)}, '
                        f'Precision: {sum(Precision) / len(Precision)}, Recall: {sum(Recall) / len(Recall)}, '
                        f'F1: {sum(F1) / len(F1)}, APscore: {sum(APscore) / len(APscore)}, Normalized_AP: {norm_ap:.4f}, '
                        f'Ami: {sum(Ami) / len(Ami)}, Ari: {sum(Ari) / len(Ari)}')

        model_dir = os.path.join(save_path, 'model_save', proj_name)
        if not os.path.exists(model_dir):
            os.makedirs(model_dir)
            print('model dir for ' + proj_name + ' has been created.')

        model_save_path = os.path.join(model_dir, 'model.ckpt')
        torch.save(model.state_dict(), model_save_path)


def set_seed(seed=42):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


if __name__ == '__main__':
    run_parser = argparse.ArgumentParser(description='Run DrugFormer on a prepared dataset')
    run_parser.add_argument('--project_name', type=str, required=True,
                            help='Dataset/project name. The dataset should be located at data_path/project_name')
    run_parser.add_argument('--data_path', type=str, required=True,
                            help='Root directory containing the prepared HuggingFace datasets')
    run_parser.add_argument('--save_path', type=str, required=True,
                            help='Directory to save trained models')
    run_parser.add_argument('--log_folder', type=str, required=True,
                            help='Directory to save log files')

    run_args, unknown = run_parser.parse_known_args()

    if not os.path.exists(run_args.log_folder):
        os.makedirs(run_args.log_folder, exist_ok=True)
        print(f'Log folder created: {run_args.log_folder}')

    if not os.path.exists(run_args.save_path):
        os.makedirs(run_args.save_path, exist_ok=True)
        print(f'Model save folder created: {run_args.save_path}')

    set_project_logger(run_args.project_name, log_folder=run_args.log_folder)
    print('Processing dataset:', str(run_args.project_name))

    set_seed(seed=42)
    run(
        proj_name=run_args.project_name,
        data_path=run_args.data_path,
        save_path=run_args.save_path
    )

    print(f'Dataset {run_args.project_name} completed.')
