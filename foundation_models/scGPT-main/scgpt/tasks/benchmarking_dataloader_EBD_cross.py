from torch.utils.data import Dataset, DataLoader
import re
import torch
import pickle

import os
import numpy as np


class BioDataset(Dataset):
    def __init__(self, data_file):
        super(BioDataset, self).__init__()
        print('loading dataset...')
        ebd_all = []
        lbl_all = []

        with open(data_file, 'rb') as pklfile:
            PKLdata = pickle.load(pklfile)

        ebd = PKLdata['embeds']
        lbl = PKLdata['labels']
        lbl_onehot = np.eye(2)[lbl]
        ebd_all.extend(torch.tensor(ebd))
        lbl_all.extend(torch.tensor(lbl_onehot))

        self.embeds = ebd_all
        self.labels = lbl_all
        assert len(ebd_all) == len(lbl_all)

        self.length = len(self.embeds)
        print('number of samples:', self.length)

    def __getitem__(self, item):
        return self.embeds[item], self.labels[item], item

    def __len__(self):
        return self.length


def dataloader(args):
    print('Loading training dataset:')
    train_dataset = BioDataset(args.train_file)

    print('Loading testing dataset:')
    test_dataset = BioDataset(args.test_file)

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

    print('train samples:', len(train_dataset))
    print('test samples:', len(test_dataset))

    return train_data_loader, test_data_loader