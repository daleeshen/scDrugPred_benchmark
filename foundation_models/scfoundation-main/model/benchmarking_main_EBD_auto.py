import torch
from benchmarking_dataloader_EBD_auto import dataloader
import argparse
from tqdm import tqdm
from torch import optim
from tool import (
    Accuracy_score,
    F1_score,
    AUROC_score,
    Recall_score,
    Precision_score
)
import pandas as pd
import pickle
from torch import nn
import csv
from pathlib import Path


class MLP_Classifier(nn.Module):
    def __init__(
        self,
        input_size=3072,
        hidden_size=512,
        num_classes=2
    ):
        super(MLP_Classifier, self).__init__()

        self.fc1 = nn.Linear(
            input_size,
            hidden_size
        )

        self.relu = nn.ReLU()

        self.fc2 = nn.Linear(
            hidden_size,
            hidden_size
        )

        self.fc3 = nn.Linear(
            hidden_size,
            num_classes
        )

    def forward(self, x):
        x = self.fc1(x)
        x = self.relu(x)

        x = self.fc2(x)
        x = self.relu(x)

        x = self.fc3(x)

        return x


###############################################################################
# Arguments
###############################################################################

parser = argparse.ArgumentParser(
    description='AI4Bio'
)

# training parameters
parser.add_argument(
    '--ep_num',
    type=int,
    default=10,
    help='epoch number of training'
)

parser.add_argument(
    '--train_batch_size',
    type=int,
    default=128,
    help=''
)

parser.add_argument(
    '--test_batch_size',
    type=int,
    default=256,
    help=''
)

parser.add_argument(
    '--data_file',
    type=str,
    required=True,
    help='Path to one *_embeds.npy file'
)

parser.add_argument(
    '--label_file',
    type=str,
    required=True,
    help='Path to the corresponding *_labels.npy file'
)

parser.add_argument(
    '--input_file',
    type=str,
    required=True,
    help='Path to the original TXT file for reading cell barcodes'
)

parser.add_argument(
    '--lr',
    type=int,
    default=0.0001,
    help=''
)

parser.add_argument(
    '--train_rate',
    type=float,
    default=0.8,
    help=''
)

parser.add_argument(
    '--ft_list',
    type=list,
    default=['none'],
    help=''
)

parser.add_argument(
    '--n_splits',
    type=int,
    default=3,
    help='Number of folds for stratified K-fold'
)

parser.add_argument(
    '--output_file',
    type=str,
    required=True,
    help='CSV file for saving final OOF metrics'
)

args = parser.parse_args()


###############################################################################
# Run
###############################################################################

def run():
    seed = 24

    torch.manual_seed(seed)

    device = (
        'cuda'
        if torch.cuda.is_available()
        else 'cpu'
    )

    ###########################################################################
    # Generate three-fold dataloaders
    ###########################################################################

    fold_loaders = dataloader(args)

    # OOF results collected from all folds
    oof_pred_all = []
    oof_lbl_all = []
    oof_ebd_all = []
    oof_idx_all = []
    oof_fold_all = []

    # Training metrics from the last epoch of each fold
    train_acc_all = []
    train_auc_all = []
    train_pre_all = []
    train_rcl_all = []
    train_f1_all = []

    ###########################################################################
    # Three-fold cross-validation
    ###########################################################################

    for fold, (
        train_data_loader,
        test_data_loader
    ) in enumerate(fold_loaders):

        print(
            '\n'
            '============================================================'
        )

        print(
            'Fold:',
            fold + 1,
            '/',
            args.n_splits
        )

        print(
            '============================================================'
        )

        #######################################################################
        # Initialize a new model for each fold
        #######################################################################

        torch.manual_seed(seed)

        model = MLP_Classifier().to(device)

        loss_function = torch.nn.BCEWithLogitsLoss()

        model_parameters = [
            parameter
            for parameter in model.parameters()
            if parameter.requires_grad
        ]

        optimizer = optim.AdamW(
            model_parameters,
            lr=args.lr
        )

        #######################################################################
        # Epochs
        #######################################################################

        for epoch in range(args.ep_num):
            ###################################################################
            # Training
            ###################################################################

            loss_sum = 0
            pred_all = []
            lbl_all = []

            with tqdm(
                train_data_loader,
                ncols=80,
                position=0,
                leave=True
            ) as batches:

                for b in batches:
                    input_embeds, labels, indices = b

                    pred = model(
                        input_embeds.float().to(device)
                    )

                    loss = loss_function(
                        pred,
                        labels.to(device)
                    )

                    optimizer.zero_grad()
                    loss.backward()
                    optimizer.step()

                    loss_sum = loss_sum + loss

                    pred_all.extend(pred)
                    lbl_all.extend(labels)

            pred_all_ = torch.stack(pred_all)
            lbl_all_ = torch.stack(lbl_all)

            acc = Accuracy_score(
                pred_all_,
                lbl_all_
            )

            f1 = F1_score(
                pred_all_,
                lbl_all_
            )

            try:
                aur = AUROC_score(
                    pred_all_,
                    lbl_all_
                )
            except:
                aur = 0.0

            pre = Precision_score(
                pred_all_,
                lbl_all_
            )

            rcl = Recall_score(
                pred_all_,
                lbl_all_
            )

            print(
                'Fold:',
                fold + 1,
                'Training epoch:',
                epoch,
                'loss:',
                loss_sum / len(batches),
                'Accuracy:',
                round(acc, 4),
                'AUROC:',
                round(aur, 4),
                'Precision:',
                round(pre, 4),
                'Recall:',
                round(rcl, 4),
                'F1:',
                round(f1, 4)
            )

            train_acc = acc
            train_auc = aur
            train_pre = pre
            train_rcl = rcl
            train_f1 = f1

            ###################################################################
            # Testing
            ###################################################################

            loss_sum = 0
            pred_all = []
            lbl_all = []
            ebd_all = []
            idx_all = []

            with torch.no_grad():
                with tqdm(
                    test_data_loader,
                    ncols=80,
                    position=0,
                    leave=True
                ) as batches:

                    for b in batches:
                        input_embeds, labels, indices = b

                        pred = model(
                            input_embeds.float().to(device)
                        )

                        loss = loss_function(
                            pred,
                            labels.to(device)
                        )

                        loss_sum = loss_sum + loss

                        pred_all.extend(pred)
                        lbl_all.extend(labels)
                        ebd_all.extend(
                            input_embeds.to('cpu')
                        )

                        idx_all.extend(
                            indices.tolist()
                        )

                pred_all_ = torch.stack(pred_all)
                lbl_all_ = torch.stack(lbl_all)
                ebd_all_ = torch.stack(ebd_all)

                acc = Accuracy_score(
                    pred_all_,
                    lbl_all_
                )

                f1 = F1_score(
                    pred_all_,
                    lbl_all_
                )

                try:
                    aur = AUROC_score(
                        pred_all_,
                        lbl_all_
                    )
                except:
                    aur = 0.0

                pre = Precision_score(
                    pred_all_,
                    lbl_all_
                )

                rcl = Recall_score(
                    pred_all_,
                    lbl_all_
                )

                print(
                    'Fold:',
                    fold + 1,
                    'Testing epoch:',
                    epoch,
                    'loss:',
                    loss_sum / len(batches),
                    'Accuracy:',
                    round(acc, 4),
                    'AUROC:',
                    round(aur, 4),
                    'Precision:',
                    round(pre, 4),
                    'Recall:',
                    round(rcl, 4),
                    'F1:',
                    round(f1, 4)
                )

            ###################################################################
            # Save final epoch results from this fold
            ###################################################################

            if epoch == args.ep_num - 1:
                train_acc_all.append(
                    float(train_acc)
                )

                train_auc_all.append(
                    float(train_auc)
                )

                train_pre_all.append(
                    float(train_pre)
                )

                train_rcl_all.append(
                    float(train_rcl)
                )

                train_f1_all.append(
                    float(train_f1)
                )

                oof_pred_all.append(
                    pred_all_.detach().cpu()
                )

                oof_lbl_all.append(
                    lbl_all_.detach().cpu()
                )

                oof_ebd_all.append(
                    ebd_all_.detach().cpu()
                )

                oof_idx_all.extend(
                    idx_all
                )

                oof_fold_all.extend(
                    [fold + 1] * len(idx_all)
                )

    ###########################################################################
    # Combine predictions from all folds
    ###########################################################################

    pred_all_ = torch.cat(
        oof_pred_all,
        dim=0
    )

    lbl_all_ = torch.cat(
        oof_lbl_all,
        dim=0
    )

    ebd_all_ = torch.cat(
        oof_ebd_all,
        dim=0
    )

    ###########################################################################
    # Calculate pooled OOF metrics
    ###########################################################################

    acc = Accuracy_score(
        pred_all_,
        lbl_all_
    )

    f1 = F1_score(
        pred_all_,
        lbl_all_
    )

    try:
        aur = AUROC_score(
            pred_all_,
            lbl_all_
        )
    except:
        aur = 0.0

    pre = Precision_score(
        pred_all_,
        lbl_all_
    )

    rcl = Recall_score(
        pred_all_,
        lbl_all_
    )

    print(
        '\n'
        '============================================================'
    )

    print('OOF results')

    print(
        '============================================================'
    )

    print(
        'Number of cells:',
        len(oof_idx_all)
    )

    print(
        'Accuracy:',
        round(acc, 4),
        'AUROC:',
        round(aur, 4),
        'Precision:',
        round(pre, 4),
        'Recall:',
        round(rcl, 4),
        'F1:',
        round(f1, 4)
    )

    print(
        '============================================================'
    )

    roc = {}

    roc['pred'] = pred_all_
    roc['label'] = lbl_all_
    roc['ebd'] = ebd_all_

    ###########################################################################
    # Check that all cells received exactly one OOF prediction
    ###########################################################################

    barcode_df = pd.read_csv(
        args.input_file,
        sep='\t',
        usecols=[0]
    )

    cell_barcode = (
        barcode_df.iloc[:, 0]
        .astype(str)
        .tolist()
    )

    if len(oof_idx_all) != len(cell_barcode):
        raise ValueError(
            'The number of OOF predictions does not match '
            'the number of cells in the original TXT file: '
            f'{len(oof_idx_all)} predictions versus '
            f'{len(cell_barcode)} cells'
        )

    if sorted(oof_idx_all) != list(
        range(len(cell_barcode))
    ):
        raise ValueError(
            'Some cells are missing or duplicated in '
            'the OOF predictions'
        )

    ###########################################################################
    # Save metrics
    ###########################################################################

    output_file = Path(
        args.output_file
    )

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    result = {
        'dataset': Path(args.data_file).stem,
        'n_splits': args.n_splits,
        'train_accuracy': (
            sum(train_acc_all)
            / len(train_acc_all)
        ),
        'train_auc': (
            sum(train_auc_all)
            / len(train_auc_all)
        ),
        'train_precision': (
            sum(train_pre_all)
            / len(train_pre_all)
        ),
        'train_recall': (
            sum(train_rcl_all)
            / len(train_rcl_all)
        ),
        'train_f1': (
            sum(train_f1_all)
            / len(train_f1_all)
        ),
        'test_accuracy': float(acc),
        'test_auc': float(aur),
        'test_precision': float(pre),
        'test_recall': float(rcl),
        'test_f1': float(f1)
    }

    pd.DataFrame(
        [result]
    ).to_csv(
        output_file,
        index=False
    )

    print(
        'Saved metrics:',
        output_file
    )

    ###########################################################################
    # Save complete cell-level OOF results table
    ###########################################################################

    predicted_score = torch.sigmoid(
        pred_all_[:, 1]
    ).detach().cpu().numpy()

    predicted_label = torch.argmax(
        pred_all_,
        dim=1
    ).detach().cpu().numpy()

    ground_truth_label = torch.argmax(
        lbl_all_,
        dim=1
    ).detach().cpu().numpy()

    prediction_result = pd.DataFrame({
        'Cell_barcode': [
            cell_barcode[i]
            for i in oof_idx_all
        ],
        'ground_truth_label': ground_truth_label,
        'prediction_score': predicted_score,
        'predicted_label': predicted_label
    })

    # Restore the original cell order in the TXT file
    prediction_result['_original_index'] = (
        oof_idx_all
    )

    prediction_result = (
        prediction_result
        .sort_values('_original_index')
        .drop(columns=['_original_index'])
        .reset_index(drop=True)
    )

    prediction_file = (
        output_file.parent
        / (
            output_file.stem
            + '_predictions.csv'
        )
    )

    prediction_result.to_csv(
        prediction_file,
        index=False
    )

    print(
        'Saved predictions:',
        prediction_file
    )

    print(
        'Number of prediction rows:',
        len(prediction_result)
    )


if __name__ == '__main__':
    run()