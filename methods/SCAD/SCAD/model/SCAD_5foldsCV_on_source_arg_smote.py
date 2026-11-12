import argparse
import torch
import numpy as np
import pandas as pd
from numpy import random
from sklearn.metrics import *
import sklearn.preprocessing as sk
import os
import sys
#from best_param_extractor import extract_best_auc_params_smote
from torch.optim.lr_scheduler import ReduceLROnPlateau
from torch.utils.data.sampler import WeightedRandomSampler
import random
from sklearn.model_selection import StratifiedKFold
import itertools
from itertools import cycle
from SCADmodules import *
####use SMOTE
from imblearn.over_sampling import SMOTE
from imblearn.pipeline import Pipeline
from imblearn.under_sampling import RandomUnderSampler
from collections import Counter
import pandas as pd
import re

##use GPU##
GPU = True

if GPU:
    device = "cuda"
else:
    device = "cpu"

init_seed = 42

#######################################################
#                    FUNCTIONS                        #
#######################################################
def folder_checker(path):
    if not os.path.exists(path):
        os.makedirs(path)
        print("Directory ", path, " Created ")
    else:
        print("Directory ", path, " already exists")

def predict_label(XTestCells, gen_model, map_model):
    """
    Inputs:
    :param XTestCells - X_target_test
    :param gen_model    - current FX model
    :param map_model    - current MTL model

    Output:
        - Predicted (binary) labels (Y for input target test data)
    """

    gen_model.eval()
    gen_model.to(device)
    map_model.eval()
    map_model.to(device)

    if extractor == 'FX' or extractor == 'FX_MLP':
        F_xt_test = gen_model(XTestCells.to(device))

    else:
        raise NotImplementedError('{} is not a valid extractor'.format(extractor))

    yhatt_test = map_model(F_xt_test)
    return yhatt_test

def evaluate_model(XTestCells, YTestCells, gen_model, map_model):
    """
    Inputs:
    :param XTestCells - Single cell test data
    :param YTestCells - true class labels (binary) for Single cell test data
    :param path_to_models - path to the saved models from training

    Outputs:
        - test loss
        - test accuracy (AUC)
    """
    XTestCells = XTestCells.to(device)
    YTestCells = YTestCells.to(device)
    y_predicted = predict_label(XTestCells, gen_model, map_model)

    # #LOSSES
    C_loss_eval = torch.nn.BCELoss()
    closs_test = C_loss_eval(y_predicted, YTestCells)

    if device == "cuda":
        YTestCells = YTestCells.to("cpu")
    yt_true_test = YTestCells.view(-1, 1)
    yt_true_test = yt_true_test.cpu()
    y_predicted = y_predicted.cpu()

    AUC_test = roc_auc_score(yt_true_test.detach().numpy(), y_predicted.detach().numpy())
    return closs_test, AUC_test

def roc_auc_score_trainval(y_true, y_predicted):
    # To handle the case where we only have training samples of one class
    # in our mini-batch when training since roc_auc_score 'breaks' when
    # there is only one class present in y_true:
    # ValueError: Only one class present in y_true. ROC AUC score is not defined in that case.

    # The following code is taken from
    # https://stackoverflow.com/questions/45139163/roc-auc-score-only-one-class-present-in-y-true?rq=1 #
    if len(np.unique(y_true)) == 1:
        return accuracy_score(y_true, np.rint(y_predicted))
    return roc_auc_score(y_true, y_predicted)

def extract_best_auc_params_smote(log_text, drug):
    # Regular expression to extract necessary information
    pattern = re.compile(
        rf"Drug = {drug}, extractor = (\w+), source average validation AUC of iter(\d+) = ([\d.]+), "
        r"hdim(\d+)_zdim(\d+)_lr([\d.]+)_epoch(\d+)_lambda([\d.]+)_"
        r"dropouts([\d.]+)_([\d.]+)_([\d.]+)_mbS(\d+)_mbT(\d+)"
    )

    best_iter = None
    best_auc = float('-inf')
    best_params = {}

    # Iterate through matches and find the best AUC
    for match in pattern.finditer(log_text):
        extractor, iter_num, auc, hdim, zdim, lr, epoch, lambda1, d1, d2, d3, mbS, mbT = match.groups()
        auc = float(auc)
        if auc > best_auc:
            best_auc = auc
            best_iter = iter_num
            best_params = {
                "Drug": drug,
                "Extractor": extractor,
                "AUC": auc,
                "hdim": int(hdim),
                "zdim": int(zdim),
                "lr": float(lr),
                "epoch": int(epoch),
                "lambda1": float(lambda1),
                "dropouts": (float(d1), float(d2), float(d3)),
                "mbS": int(mbS),
                "mbT": int(mbT),
                "iter": int(iter_num)
            }

    return best_params




# dataset name list
GSE_folder = ['GSE108397_PLX4720', 'GSE115251_NVP-TAE684', 'GSE131984_JQ1_DMSO_rm_JQ1', 'GSE131984_Pac_DMSO_rm_Paclitaxel', 'GSE131984_Pal_DMSO_rm_PD-0332991', 'GSE140440_Docetaxel', 'GSE149214_Erlotinib', 'GSE158457_Dasatinib', 'GSE163836_Paclitaxel', 'GSE169246_PacBlood_Paclitaxel', 'GSE169246_PacTissue_Paclitaxel', 'GSE175716_Sorafenib', 'GSE186960_Gemcitabine']

#遍历数据集
for folder_name in GSE_folder:
    GEO_id = folder_name.rsplit('_', 1)[0]
    Drug = folder_name.rsplit('_', 1)[1]
    DRUG = folder_name.rsplit('_', 1)[1]
    extractor = "FX"
    geneset = ['_norm','_ppi', '_tp4k']
    adjust = 'N'
    init_seed = 42
    if GEO_id == 'GSE158457':
        cell_lines = 'all'
    else:
        cell_lines = 'solid'

    random.seed(init_seed)
    np.random.seed(init_seed)
    os.environ['PYTHONHASHSEED'] = str(init_seed)
    torch.manual_seed(init_seed)
    torch.cuda.manual_seed(init_seed)
    torch.cuda.manual_seed_all(init_seed)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True

    for gs in geneset:
        feature = gs
        LOG_FILE_PATH = 'D:/SCAD/results_smote' + gs + '/train_eval_on_source_hypers_summary_5fds_smote_fixlr.txt'
        with open(LOG_FILE_PATH, 'r') as log_file:
            log_content = log_file.read()
        best_params = extract_best_auc_params_smote(log_content, drug=DRUG)

        h_dim = best_params['hdim']
        z_dim = best_params['zdim']
        epoch = best_params['epoch']
        lam1 = best_params['lambda1']
        mbS = best_params['mbS']
        mbT = best_params['mbT']

        if adjust == "Y":
            X_5k_cells_raw = pd.read_csv("D:/SCAD/data/Pancan_CCLE_10X/scrna_ccle_CPM_random_5k_raw_x_scaled_ncbi.tsv",
                                         sep='\t', index_col=0, decimal='.')

            feature = geneset + '_adj'

        tr_ts_summary_path = "D:/SCAD/best_param_results_smote" + feature + "/" + GEO_id + "_" + DRUG + "_" + "train_test_summary_5folds_5seeds.txt"
        max_iter = 1
        MODE = "bimodal/best_models/" + extractor
        SAVE_RESULTS_TO = "D:/SCAD/best_param_results_smote" + feature + "/" + folder_name + "/" + MODE + "/"
        SAVE_TRACE_TO = "D:/SCAD/best_param_results_smote" + feature + "/" + folder_name + "/" + MODE + "/trace/"

        if gs == '_norm':
            SOURCE_DIR = 'source_5_folds' + gs
            TARGET_DIR = 'target_5_folds' + gs
        else:
            SOURCE_DIR = 'source_5_folds_' + gs
            TARGET_DIR = 'target_5_folds_' + gs

        LOAD_DATA_FROM = 'D:/SCAD/data/split_norm/' + folder_name + '/stratified'+ gs +'/'  #####Load data from
        dirName = SAVE_RESULTS_TO + 'test/'
        folder_checker(dirName)

        #######################################################
        #                Hyper-Parameter Lists                #
        #######################################################

        ls_splits = ['split1', 'split2', 'split3', 'split4', 'split5']
        ls_mb_size = [{'mbS': mbS, 'mbT': mbT}]  # Set the batch size for the best model manually

        if gs == '_norm':
            YGDSC = pd.read_csv("D:/SCAD/data/split_norm/" + folder_name + '/' + 'Source_exprs_resp_z.' + DRUG + '.tsv',
                                sep='\t', index_col=0, decimal='.')
        else:
            YGDSC = pd.read_csv("D:/SCAD/data/split_norm/" + folder_name + '/' + 'Source_exprs_resp_z.' + DRUG + gs + '.tsv',
                                sep='\t', index_col=0, decimal='.')

        threshold = YGDSC['logIC50'][YGDSC['response'] == 0].min()
        Counter(YGDSC['response'])[0] / len(YGDSC['response'])
        Counter(YGDSC['response'])[1] / len(YGDSC['response'])
        class_sample_count = np.array([Counter(YGDSC['response'])[0] / len(YGDSC['response']),
                                       Counter(YGDSC['response'])[1] / len(YGDSC['response'])])

        #######################################################
        #         AITL Model Training Starts Here             #
        #######################################################
        for index, mbsize in enumerate(ls_mb_size):
            mbS = mbsize['mbS']
            mbT = mbsize['mbT']

            # Hard-code the values for testing and training 'best model' ##
            lr = 1e-3
            lambdvae = 0
            dropout_gen = 0.5
            dropout_mtl = dropout_gen
            dropout_dg = dropout_gen

            print("-- Parameters used: --")
            print("h_dim: {}\nz_dim: {}\nlr: {}\nepoch: {}\nlambda1: {}\n".format(h_dim, z_dim, lr, epoch, lam1))

            print("mbS: {}\nmbT: {}\ndropout_gen: {}\ndropout_mtl: {}\ndropout_dg: {}\n".format(mbS,
                                                                                                mbT,
                                                                                                dropout_gen,
                                                                                                dropout_mtl,
                                                                                                dropout_dg))

            batch_sizes = 'mbS' + str(mbS) + '_mbT' + str(mbT)
            test_results_name = 'hdim' + str(h_dim) + '_zdim' + str(z_dim) + '_lr' + str(lr) + '_epoch' + str(
                epoch) + '_lambda1' + str(lam1) + \
                                '_mbS' + str(mbS) + '_mbT' + str(mbT) + '_smote' + '_seed' + str(init_seed) + '.tsv'
            test_predicted_y_name = 'hdim' + str(h_dim) + '_zdim' + str(z_dim) + '_lr' + str(lr) + '_epoch' + str(
                epoch) + '_lambda1' + str(lam1) + \
                                    '_mbS' + str(mbS) + '_mbT' + str(mbT) + '_seed' + str(
                init_seed) + '_predicted_y.tsv'


            test_results_dir = SAVE_RESULTS_TO + 'test/' + batch_sizes + '/'
            dirName = test_results_dir
            folder_checker(dirName)

            test_results_file = os.path.join(test_results_dir, test_results_name)
            if os.path.isfile(test_results_file):
                os.remove(test_results_file)

            test_predicted_y_file = os.path.join(test_results_dir, test_predicted_y_name)
            if os.path.isfile(test_predicted_y_file):
                os.remove(test_predicted_y_file)

            with open(test_results_file, 'a') as f:
                f.write("-- Parameters --\n\n")
                f.write(
                    "h_dim: {}\nz_dim: {}\nlr: {}\nepoch: {}\nlambda1: {}\nmbS: {}\nmbT: {}\n".format(h_dim, z_dim, lr,
                                                                                                      epoch,
                                                                                                      lam1, mbS, mbT))

            AUCtest_splits_total = []
            APRtest_splits_total = []

            for split in ls_splits:  # for each split
                print("\n\nReading data for {} ...\n".format(split))
                # Loading Source Data, in this model, both source and target label are binarized, the target domain labels are masked during training#
                XTrainGDSC = pd.read_csv(LOAD_DATA_FROM + SOURCE_DIR + '/' + split + '/X_train_source.tsv',
                                         sep='\t', index_col=0, decimal='.')
                YTrainGDSC = pd.read_csv(LOAD_DATA_FROM + SOURCE_DIR + '/' + split + '/Y_train_source.tsv',
                                         sep='\t', index_col=0, decimal='.')
                XValGDSC = pd.read_csv(LOAD_DATA_FROM + SOURCE_DIR + '/' + split + '/X_val_source.tsv',
                                       sep='\t', index_col=0, decimal='.')
                YValGDSC = pd.read_csv(LOAD_DATA_FROM + SOURCE_DIR + '/' + split + '/Y_val_source.tsv',
                                       sep='\t', index_col=0, decimal='.')

                try:
                    over = SMOTE(sampling_strategy=0.5)
                    under = RandomUnderSampler(sampling_strategy=0.5)
                    steps = [('o', over), ('u', under)]
                    pipeline = Pipeline(steps=steps)
                    ##
                    ##conduct SMOTE oversampling on source training set
                    XTrainGDSC, YTrainGDSC = pipeline.fit_resample(XTrainGDSC, YTrainGDSC)
                except ValueError:
                    over = SMOTE(sampling_strategy=0.6)
                    under = RandomUnderSampler(sampling_strategy=0.6)
                    steps = [('o', over), ('u', under)]
                    pipeline = Pipeline(steps=steps)
                    ##
                    ##conduct SMOTE oversampling on source training set
                    XTrainGDSC, YTrainGDSC = pipeline.fit_resample(XTrainGDSC, YTrainGDSC)

                # Loading Target (Single cell) Data, YTrainCells are not used during training#
                XTrainCells = pd.read_csv(LOAD_DATA_FROM + TARGET_DIR + '/' + split + '/X_train_target.tsv',
                                          sep='\t', index_col=0, decimal='.')
                YTrainCells = pd.read_csv(LOAD_DATA_FROM + TARGET_DIR + '/' + split + '/Y_train_target.tsv',
                                          sep='\t', index_col=0, decimal='.')
                XTestCells = pd.read_csv(LOAD_DATA_FROM + TARGET_DIR + '/' + split + '/X_test_target.tsv',
                                         sep='\t', index_col=0, decimal='.')
                YTestCells = pd.read_csv(LOAD_DATA_FROM + TARGET_DIR + '/' + split + '/Y_test_target.tsv',
                                         sep='\t', index_col=0, decimal='.')

                # Load randomly sampled Data #
                if adjust == "Y":
                    X_5k_cells = X_5k_cells_raw.transpose()
                    X_5k_cells.columns = X_5k_cells.columns.astype(str)
                    intersect_gene = [value for value in XValGDSC.columns if value in X_5k_cells.columns]
                    X_5k = X_5k_cells[intersect_gene]
                    XTrainGDSC = XTrainGDSC[intersect_gene]
                    XValGDSC = XValGDSC[intersect_gene]
                    XTrainCells = XTrainCells[intersect_gene]
                    XTestCells = XTestCells[intersect_gene]

                print("Data successfully read!")

                ####
                model_params = 'hdim' + str(h_dim) + '_zdim' + str(z_dim) + '_lr' + str(lr) + '_epoch' + str(
                    epoch) + '_lamb1' + str(lam1) \
                               + '_mbS' + str(mbS) + '_mbT' + str(mbT)

                ###
                # Temporarily combine Source training data and Target training data
                # to fit standard scaler on gene expression of combined training data.
                # Then, apply fitted scaler to (and transform) Source validation,
                # Target validation, and Target test (e.g. normalize validation and test
                # data of source and target with respect to source and target train)
                XTrainCombined = pd.concat([XTrainGDSC, XTrainCells])

                if adjust == "Y":
                    XTrainCombined = pd.concat([XTrainCombined, X_5k_cells[intersect_gene]])

                scalerTrain = sk.StandardScaler()
                scalerTrain.fit(XTrainCombined.values)
                # N means normalize
                XTrainGDSC_N = scalerTrain.transform(XTrainGDSC.values)
                XTrainCells_N = scalerTrain.transform(XTrainCells.values)

                XTestCells_N = scalerTrain.transform(XTestCells.values)  # scalerTrain.fit(XTrainCombined.values)
                TXTestCells_N = torch.FloatTensor(XTestCells_N)
                TYTestCells = torch.FloatTensor(YTestCells.values.astype(int))
                TXTestCells_N = TXTestCells_N.to(device)
                TYTestCells = TYTestCells.to(device)

                weight = 1. / class_sample_count
                # upsampling source domain training set which is unbalanced##
                samples_weight = np.array([weight[t] for t in YTrainGDSC.values])

                samples_weight = torch.from_numpy(samples_weight)
                samples_weight = samples_weight.reshape(-1)  # Flatten out the weights so it's a 1-D tensor of weights
                sampler = WeightedRandomSampler(samples_weight.type('torch.DoubleTensor'), len(samples_weight),
                                                replacement=True)

                # Apply sampler on XTrainCells_N
                ## GDSC cell line dataset ##
                CDataset = torch.utils.data.TensorDataset(torch.FloatTensor(XTrainGDSC_N),
                                                          torch.FloatTensor(YTrainGDSC.values))

                PDataset = torch.utils.data.TensorDataset(torch.FloatTensor(XTrainCells_N),
                                                          torch.FloatTensor(YTrainCells.values.astype(int)))
                # upsampling source domain training set which is unbalanced##
                # CLoader = torch.utils.data.DataLoader(dataset = CDataset, batch_size= mbS, shuffle=False, sampler = sampler)
                CLoader = torch.utils.data.DataLoader(dataset=CDataset, batch_size=mbS, shuffle=True)
                PLoader = torch.utils.data.DataLoader(dataset=PDataset, batch_size=mbT, shuffle=True)

                ##
                n_sample, IE_dim = XTrainGDSC_N.shape

                if extractor == 'FX':
                    Gen = FX(dropout_gen, IE_dim, h_dim, z_dim)
                else:
                    raise NotImplementedError('{} is not a valid extractor'.format(extractor))

                Map = MTLP(dropout_mtl, h_dim, z_dim)  #### Mapper
                DG = Discriminator(dropout_dg, h_dim, z_dim)  #### Global Discriminator####

                Gen.to(device)
                Map.to(device)
                DG.to(device)

                # feature extractor, domain adaptation and global discriminator #
                optimizer_2 = torch.optim.Adagrad(itertools.chain(Gen.parameters(), Map.parameters(), DG.parameters()),
                                                  lr=lr)

                C_loss = torch.nn.BCELoss()

                l1 = []
                l2 = []
                classif = []
                L = []  # total loss
                DG_losstr = []
                DG_auctr = []

                for it in range(epoch):
                    epoch_cost1ls = []
                    epoch_cost2ls = []
                    epoch_auc = []
                    epoch_classifls = []
                    epoch_DGloss = []
                    epoch_DGauc = []

                    epoch_loss = []

                    for i, data in enumerate(zip(cycle(CLoader), PLoader)):
                        DataS = data[0]
                        DataT = data[1]
                        #
                        xs = DataS[0].to(device)
                        ys = DataS[1].view(-1, 1).to(device)
                        xt = DataT[0].to(device)
                        yt = DataT[1].view(-1, 1).to(device)

                        # Skip to next set of training batch if any of xs or xt has less
                        # than a certain threshold of training examples.
                        # Let such threshold be 5
                        if xs.size()[0] < 5 or xt.size()[0] < 5:
                            continue

                        Gen.train()  ### Gen = EN(dropout_gen, IE_dim, h_dim, z_dim)
                        Map.train()  ### mapper
                        DG.train()  ### discriminator

                        if extractor == 'FX' or extractor == 'FX_SFMX' or extractor == 'FX_MLP':
                            F_xs = Gen(xs)
                            F_xt = Gen(xt)

                        else:
                            raise NotImplementedError('{} is not a valid extractor'.format(extractor))

                        yhat_xs = Map(F_xs)
                        # classification for GDSC
                        closs = C_loss(yhat_xs, ys)
                        loss1 = closs

                        Labels = torch.ones(F_xs.size(0), 1)
                        Labelt = torch.zeros(F_xt.size(0), 1)
                        Lst = torch.cat([Labels, Labelt], 0).to(device)  ## combine domain labels
                        Xst = torch.cat([F_xs, F_xt], 0).to(device)  ##combine domain data

                        yhat_DG = DG(Xst)  ##predicted domain label from global discriminator
                        DG_loss = C_loss(yhat_DG, Lst)

                        loss2 = lam1 * DG_loss

                        Loss = loss1 + loss2

                        ###

                        optimizer_2.zero_grad()
                        Loss.backward()
                        optimizer_2.step()

                        epoch_cost1ls.append(loss1)
                        epoch_cost2ls.append(loss2)
                        epoch_classifls.append(closs)
                        epoch_loss.append(Loss)
                        epoch_DGloss.append(DG_loss)

                        y_true = yt.view(-1, 1)
                        y_true = y_true.cpu()

                        y_trueDG = Lst.view(-1, 1)
                        y_predDG = yhat_DG

                        y_trueDG = y_trueDG.cpu()
                        y_predDG = y_predDG.cpu()

                        AUCDG = roc_auc_score_trainval(y_trueDG.detach().numpy(), y_predDG.detach().numpy())
                        epoch_DGauc.append(AUCDG)

                    l1.append(torch.mean(torch.Tensor(epoch_cost1ls)))
                    l2.append(torch.mean(torch.Tensor(epoch_cost2ls)))
                    classif.append(torch.mean(torch.Tensor(epoch_classifls)))
                    L.append(torch.mean(torch.FloatTensor(epoch_loss)))
                    DG_losstr.append(torch.mean(torch.Tensor(epoch_DGloss)))
                    DG_auctr.append(torch.mean(torch.Tensor(epoch_DGauc)))

                    # Take average across all training batches
                    loss1tr_mean = l1[it]
                    loss2tr_mean = l2[it]
                    totlosstr_mean = L[it]
                    DGlosstr_mean = DG_losstr[it]
                    DGauctr_mean = DG_auctr[it]
                    closstr_mean = classif[it]

                    print("\n\nEpoch: {}".format(it))
                    print("(train) loss1 mean: {}".format(loss1tr_mean))
                    print("(train) loss2 mean: {}".format(loss2tr_mean))
                    print("(train) total loss mean: {}".format(totlosstr_mean))
                    print("################totlosstr_mean={}###############".format(totlosstr_mean))
                    print("(train) DG loss mean: {}".format(DGlosstr_mean))

                    # Write to file
                    # Take avg
                    save_model_to = SAVE_TRACE_TO + batch_sizes + '/' + model_params + '/model/'
                    if not os.path.exists(save_model_to):
                        os.makedirs(save_model_to)
                        print("Directory ", save_model_to, " Created ")

                    save_best_model_to = os.path.join(save_model_to, split + '_best_model.pt')
                    ## Saving multiple models in one file (Feature extractors, mapper, discriminators) ##
                    ##  https://pytorch.org/tutorials/beginner/saving_loading_models.html#saving-multiple-models-in-one-file ##
                    # torch.save({
                    #            'Gen_state_dict': Gen.state_dict(),
                    #            'Map_state_dict': Map.state_dict(),
                    #            'DG_state_dict': DG.state_dict(),
                    #            'optimizer_2_state_dict': optimizer_2.state_dict(),
                    #            }, save_best_model_to)

                ## Evaluate model ##
                print("\n\n-- Evaluation -- \n\n")
                print("TXTestCells_N shape: {}\n".format(TXTestCells_N.size()))
                print("TYTestCells shape: {}\n".format(TYTestCells.size()))
                save_best_model_to = os.path.join(save_model_to, split + '_best_model.pt')
                test_loss, test_auc, test_apr, test_predicted_y = evaluate_model(TXTestCells_N, TYTestCells, Gen,
                                                                                 Map)

                print("test_apr = {}\n".format(test_apr))
                print("\n\n-- Test Results --\n\n")
                print("test loss: {}".format(test_loss))
                print("test auc: {}".format(test_auc))
                print("\n ----------------- \n\n\n")
                with open(test_results_file, 'a') as f:
                    f.write("-- Split {} --\n".format(split))
                    f.write("Test loss: {}\t Test AUC: {}\t Test APR: {}\n\n\n".format(test_loss, test_auc, test_apr))
                AUCtest_splits_total.append(test_auc)
                APRtest_splits_total.append(test_apr)

                test_predicted_y_df = pd.DataFrame(test_predicted_y.cpu().detach().numpy(), index=YTestCells.index,
                                                   columns=['Predicted_Y'])

                test_predicted_y_file = os.path.join(test_results_dir, test_predicted_y_name)
                test_predicted_y_df.to_csv(test_predicted_y_file, sep='\t')

                print(f"Predicted Y saved to {test_predicted_y_file}")

            ## Calculate Test set's avg AUC across different splits
            AUCtest_splits_total = np.array(AUCtest_splits_total)
            APRtest_splits_total = np.array(APRtest_splits_total)
            avgAUC = np.mean(AUCtest_splits_total)
            stdAUC = np.std(AUCtest_splits_total)
            avgAPR = np.mean(APRtest_splits_total)
            stdAPR = np.std(APRtest_splits_total)

            ####
            with open(test_results_file, 'a') as f:
                f.write("\n\n-- Average Test AUC --\n\n")
                f.write("Mean: {}\tStandard Deviation: {}\n".format(avgAUC, stdAUC))
                f.write("\n-- Average Test Precision Recall --\n\n")
                f.write("Mean: {}\tStandard Deviation: {}\n".format(avgAPR, stdAPR))

            with open(tr_ts_summary_path, "a") as f:
                f.write("Average Test AUC = {}; Average Test Precision Recall = {};  path = {}\n".format(avgAUC, avgAPR,
                                                                                                         test_results_file))
                f.close()
                ##


