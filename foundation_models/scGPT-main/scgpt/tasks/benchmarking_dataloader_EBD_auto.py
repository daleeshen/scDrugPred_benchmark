from torch.utils.data import Dataset, DataLoader, Subset
from sklearn.model_selection import StratifiedKFold
import re
import torch
import pickle

import os
import numpy as np


class BioDataset(Dataset):
    def __init__(self, args):
        super(BioDataset, self).__init__()
        print('loading dataset...')
        ebd_all = []
        lbl_all = []

        # Modified: load one specified pkl file instead of all files in data_path
        with open(args.data_file, 'rb') as pklfile:
            PKLdata = pickle.load(pklfile)

        ebd = PKLdata['embeds']
        lbl = PKLdata['labels']
        lbl_onehot = np.eye(2)[lbl]
        ebd_all.extend(torch.tensor(ebd))
        lbl_all.extend(torch.tensor(lbl_onehot))

        self.embeds = ebd_all
        self.labels = lbl_all
        self.raw_labels = np.array(lbl)

        assert len(ebd_all) == len(lbl_all)

        self.length = len(self.embeds)
        print('number of samples:', self.length)

    def __getitem__(self, item):
        return self.embeds[item], self.labels[item], item

    def __len__(self):
        return self.length


def dataloader(args):
    ebddataset = BioDataset(args)

    n_splits = getattr(args, 'n_splits', 3)

    skf = StratifiedKFold(
        n_splits=n_splits,
        shuffle=True,
        random_state=24
    )

    fold_loaders = []

    for fold, (train_index, test_index) in enumerate(
            skf.split(np.zeros(len(ebddataset)), ebddataset.raw_labels)):

        train_dataset = Subset(ebddataset, train_index)
        test_dataset = Subset(ebddataset, test_index)

        train_data_loader = DataLoader(
            dataset=train_dataset,
            batch_size=args.train_batch_size,
            shuffle=True
        )

        test_data_loader = DataLoader(
            dataset=test_dataset,
            batch_size=args.test_batch_size,
            shuffle=True
        )

        fold_loaders.append((train_data_loader, test_data_loader))

        print('Fold:', fold + 1,
              'train samples:', len(train_dataset),
              'test samples:', len(test_dataset))

    return fold_loaders