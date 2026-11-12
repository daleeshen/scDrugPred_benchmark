import torch
from benchmarking_dataloader_EBD_cross import dataloader
import argparse
from tqdm import tqdm
from torch import optim
from tool import Accuracy_score, F1_score, AUROC_score, Recall_score, Precision_score
import pandas as pd
import pickle
from torch import nn
import csv
from pathlib import Path


parser = argparse.ArgumentParser(description='AI4Bio')

# training parameters
parser.add_argument('--ep_num', type=int, default=10, help='epoch number of training')
parser.add_argument('--train_batch_size', type=int, default=128, help='')
parser.add_argument('--test_batch_size', type=int, default=256, help='')
parser.add_argument('--train_file', type=str, required=True, help='Path to training *_samples.pkl file')
parser.add_argument('--test_file', type=str, required=True, help='Path to testing *_samples.pkl file')
parser.add_argument('--lr', type=float, default=0.0001, help='')
parser.add_argument('--ft_list', type=list, default=['none'], help='')
parser.add_argument('--output_file', type=str, default=None, help='CSV file for saving final metrics')

args = parser.parse_args()


class MLP_Classifier(nn.Module):
    def __init__(self, input_size=512, hidden_size=512, num_classes=2):
        super(MLP_Classifier, self).__init__()
        self.fc1 = nn.Linear(input_size, hidden_size)
        self.relu = nn.ReLU()
        self.fc2 = nn.Linear(hidden_size, hidden_size)
        self.fc3 = nn.Linear(hidden_size, num_classes)

    def forward(self, x):
        x = self.fc1(x)
        x = self.relu(x)
        x = self.fc2(x)
        x = self.relu(x)
        x = self.fc3(x)
        return x


def run():
    seed = 24
    torch.manual_seed(seed)
    device = 'cuda' if torch.cuda.is_available() else 'cpu'

    model = MLP_Classifier().to(device)

    loss_function = torch.nn.BCEWithLogitsLoss()
    model_parameters = [parameter for parameter in model.parameters() if parameter.requires_grad]
    optimizer = optim.AdamW(model_parameters, lr=args.lr)

    train_data_loader, test_data_loader = dataloader(args)

    for epoch in range(args.ep_num):
        loss_sum = 0
        pred_all = []
        lbl_all = []

        with tqdm(train_data_loader, ncols=80, position=0, leave=True) as batches:
            for b in batches:
                input_embeds, labels, indices = b
                pred = model(input_embeds.float().to(device))
                loss = loss_function(pred, labels.to(device))

                optimizer.zero_grad()
                loss.backward()
                optimizer.step()

                loss_sum = loss_sum + loss
                pred_all.extend(pred)
                lbl_all.extend(labels)

        pred_all_ = torch.stack(pred_all)
        lbl_all_ = torch.stack(lbl_all)

        acc = Accuracy_score(pred_all_, lbl_all_)
        f1 = F1_score(pred_all_, lbl_all_)
        try:
            aur = AUROC_score(pred_all_, lbl_all_)
        except:
            aur = 0.0
        pre = Precision_score(pred_all_, lbl_all_)
        rcl = Recall_score(pred_all_, lbl_all_)

        print('Training epoch:', epoch, 'loss:', loss_sum/len(batches),
              'Accuracy:', round(acc, 4), 'AUROC:', round(aur, 4),
              'Precision:', round(pre, 4), 'Recall:', round(rcl, 4), 'F1:', round(f1, 4))

        train_acc = acc
        train_auc = aur
        train_pre = pre
        train_rcl = rcl
        train_f1 = f1

        loss_sum = 0
        pred_all = []
        lbl_all = []
        ebd_all = []
        idx_all = []

        with torch.no_grad():
            with tqdm(test_data_loader, ncols=80, position=0, leave=True) as batches:
                for b in batches:
                    input_embeds, labels, indices = b
                    pred = model(input_embeds.float().to(device))
                    loss = loss_function(pred, labels.to(device))

                    loss_sum = loss_sum + loss
                    pred_all.extend(pred)
                    lbl_all.extend(labels)
                    ebd_all.extend(input_embeds.to('cpu'))
                    idx_all.extend(indices.tolist())

            pred_all_ = torch.stack(pred_all)
            lbl_all_ = torch.stack(lbl_all)
            ebd_all_ = torch.stack(ebd_all)

            acc = Accuracy_score(pred_all_, lbl_all_)
            f1 = F1_score(pred_all_, lbl_all_)
            try:
                aur = AUROC_score(pred_all_, lbl_all_)
            except:
                aur = 0.0
            pre = Precision_score(pred_all_, lbl_all_)
            rcl = Recall_score(pred_all_, lbl_all_)

            print('Testing epoch:', epoch, 'loss:', loss_sum/len(batches),
                  'Accuracy:', round(acc, 4), 'AUROC:', round(aur, 4),
                  'Precision:', round(pre, 4), 'Recall:', round(rcl, 4), 'F1:', round(f1, 4))

        if epoch == args.ep_num - 1:
            roc = {}
            roc['pred'] = pred_all_
            roc['label'] = lbl_all_
            roc['ebd'] = ebd_all_

            if args.output_file is not None:
                output_file = Path(args.output_file)
                output_file.parent.mkdir(parents=True, exist_ok=True)

                result = {
                    'train_dataset': Path(args.train_file).stem,
                    'test_dataset': Path(args.test_file).stem,
                    'train_accuracy': float(train_acc),
                    'train_auc': float(train_auc),
                    'train_precision': float(train_pre),
                    'train_recall': float(train_rcl),
                    'train_f1': float(train_f1),
                    'test_accuracy': float(acc),
                    'test_auc': float(aur),
                    'test_precision': float(pre),
                    'test_recall': float(rcl),
                    'test_f1': float(f1)
                }

                pd.DataFrame([result]).to_csv(output_file, index=False)
                print('Saved metrics:', output_file)

                # Save cell-level predictions for the external test dataset
                with open(args.test_file, 'rb') as pklfile:
                    PKLdata = pickle.load(pklfile)

                cell_barcode = PKLdata['cell_barcode']
                condition = PKLdata['condition']
                predicted_score = torch.sigmoid(pred_all_[:, 1]).detach().cpu().numpy()
                predicted_label = torch.argmax(pred_all_, dim=1).detach().cpu().numpy()

                prediction_result = pd.DataFrame({
                    'Cell_barcode': [cell_barcode[i] for i in idx_all],
                    'Condition': [condition[i] for i in idx_all],
                    'predicted_score': predicted_score,
                    'predicted_label': predicted_label
                })

                prediction_file = output_file.parent / (output_file.stem + '_predictions.csv')
                prediction_result.to_csv(prediction_file, index=False)
                print('Saved predictions:', prediction_file)

            # with open('./output/testset_records.pkl', 'wb') as rocpkl:
            #     pickle.dump(roc, rocpkl)
            # torch.save(model.state_dict(), './output/model.pth')


if __name__ == '__main__':
    run()