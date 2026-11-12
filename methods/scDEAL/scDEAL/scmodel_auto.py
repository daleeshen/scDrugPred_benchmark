#!/usr/bin/env python
# coding: utf-8
import argparse
import pandas as pd
from pandas.core.frame import DataFrame
import logging
import os
import sys
import time
import numpy as np
import pandas as pd
import scanpy as sc
import torch
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from torch import nn, optim
from torch.optim import lr_scheduler
from torch.utils.data import DataLoader, TensorDataset
import DaNN.mmd as mmd
import scanpypip.preprocessing as pp
import trainers as t
import utils as ut
# from captum.attr import IntegratedGradients
from models import (AEBase, DaNN, PretrainedPredictor,
                    PretrainedVAEPredictor, VAEBase)
from scipy.spatial import distance_matrix, minkowski_distance, distance
import random
from sklearn.preprocessing import LabelEncoder


def set_seed(seed=42):
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    random.seed(seed)
    np.random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


torch.set_num_threads(20)

os.environ["OMP_NUM_THREADS"] = "20"
os.environ["MKL_NUM_THREADS"] = "20"

import warnings

warnings.filterwarnings('ignore')



class TargetModel(nn.Module):
    def __init__(self, source_predcitor, target_encoder):
        super(TargetModel, self).__init__()
        self.source_predcitor = source_predcitor
        self.target_encoder = target_encoder

    def forward(self, X_target, C_target=None):

        if (type(C_target) == type(None)):
            x_tar = self.target_encoder.encode(X_target)
        else:
            x_tar = self.target_encoder.encode(X_target, C_target)
        y_src = self.source_predcitor.predictor(x_tar)
        return y_src


def run_main(args):
    ################################################# START SECTION OF LOADING PARAMETERS #################################################
    # Read parameters
    # t0 = time.time()
    set_seed()
    # Overwrite params if checkpoint is provided
    # args.checkpoint = "save/sc_pre/integrate_data_GSE112274_drug_GEFITINIB_bottle_256_edim_512,256_pdim_256,128_model_DAE_dropout_0.1_gene_F_lr_0.5_mod_new_sam_no_DaNN.pkl"
    raw_data_path = args.raw_data_path
    root_path = args.root_path
    if (args.checkpoint not in ["False", "True"]):
        selected_model = args.checkpoint
        split_name = selected_model.split("/")[-1].split("_")
        para_names = (split_name[1::2])
        paras = (split_name[0::2])
        args.bulk_h_dims = paras[4]
        args.sc_h_dims = paras[4]
        args.predictor_h_dims = paras[5]
        args.bottleneck = int(paras[3])
        args.drug = paras[2]
        args.dropout = float(paras[7])
        args.dimreduce = paras[6]

        if (paras[0].find("GSE117872") >= 0):
            args.sc_data = "GSE117872"
            args.batch_id = paras[1].split("GSE117872")[1]
        elif (paras[0].find("MIX-Seq") >= 0):
            args.sc_data = "MIX-Seq"
            args.batch_id = paras[1].split("MIX-Seq")[1]
        else:
            args.sc_data = paras[1]

    # Laod parameters from args
    epochs = args.epochs
    dim_au_out = args.bottleneck  # 8, 16, 32, 64, 128, 256,512
    na = args.missing_value
    # if args.sc_data=='GSE117872_HN120':
    #     data_path = DATA_MAP['GSE117872']
    # elif args.sc_data=='GSE117872_HN137':
    #     data_path = DATA_MAP['GSE117872']

    data_path = args.sc_expression_file
    print("data_path =", data_path, flush=True)

    test_size = args.test_size
    if args.drug.upper() == 'PLX-4720':
        select_drug = 'PLX.4720'
    else:
        select_drug = args.drug.upper()
    freeze = args.freeze_pretrain
    valid_size = args.valid_size
    g_disperson = args.var_genes_disp
    min_n_genes = args.min_n_genes
    max_n_genes = args.max_n_genes
    log_path = os.path.join(root_path, args.logging_file)
    batch_size = args.batch_size
    encoder_hdims = args.bulk_h_dims.split(",")
    encoder_hdims = list(map(int, encoder_hdims))
    data_name = args.sc_data
    label_path = args.label
    reduce_model = args.dimreduce
    predict_hdims = args.predictor_h_dims.split(",")
    predict_hdims = list(map(int, predict_hdims))
    leiden_res = args.cluster_res
    bulk_model_path = os.path.join(root_path, args.bulk_model_path)
    sc_model_path = os.path.join(root_path, args.sc_model_path)
    sc_encoder_path = os.path.join(root_path, args.sc_encoder_path)

    # load_model = bool(args.load_sc_model)
    mod = args.mod

    # Merge parameters as string for saving model and logging
    para = str(args.bulk) + "_data_" + str(args.sc_data) + "_drug_" + str(args.drug) + "_bottle_" + str(
        args.bottleneck) + "_edim_" + str(args.bulk_h_dims) + "_pdim_" + str(
        args.predictor_h_dims) + "_model_" + reduce_model + "_dropout_" + str(args.dropout) + "_gene_" + str(
        args.printgene) + "_lr_" + str(args.lr) + "_mod_" + str(args.mod) + "_sam_" + str(args.sampling)
    source_data_path = args.bulk_data

    # Record time
    # now=time.strftime("%Y-%m-%d-%H-%M-%S")
    # Initialize logging and std out
    # out_path = log_path+now+"transfer.err"
    # log_path = log_path+now+"transfer.log"
    # out=open(out_path,"w")
    # sys.stderr=out

    # Logging parameters
    # logging.basicConfig(level=logging.INFO,
    #                 filename=log_path,
    #                 filemode='a',
    #                 format=
    #                 '%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s'
    #                 )
    # logging.getLogger('matplotlib.font_manager').disabled = True
    # logging.info(args)
    # logging.info("Start at " + str(t0))

    # Create directories if they do not exist
    for path in [log_path, bulk_model_path, sc_model_path, sc_encoder_path]:
        if not os.path.exists(path):
            # Create a new directory because it does not exist
            os.makedirs(path)
            print("The new directory is created!")

    # Save arguments
    # Overwrite params if checkpoint is provided
    if (args.checkpoint not in ["False", "True"]):
        para = os.path.basename(selected_model).split("_DaNN.pkl")[0]
        args.checkpoint = "True"

    source_model_path = args.bulk_model_path + para
    print(source_model_path)
    sc_encoder_path = args.sc_encoder_path + para
    target_model_path = args.sc_model_path + para
    # args_df = ut.save_arguments(args,now)
    ################################################# END SECTION OF LOADING PARAMETERS ##############################################################

    ################################################# START SECTION OF SINGLE CELL DATA REPROCESSING #################################################
    # Load data and preprocessing

    adata = pp.read_sc_file(data_path, meta_path=args.meta_file)

    print(f"[Check] After loading raw data: cells = {adata.n_obs}, genes = {adata.n_vars}")
    # print("sc_data", adata)
    # if data_name == 'GSE117872_HN137':
    #     adata =  ut.specific_process(adata,dataname='GSE117872',select_origin='HN137')
    # elif data_name == 'GSE117872_HN120':
    #     adata =  ut.specific_process(adata,dataname='GSE117872',select_origin='HN120')
    # elif data_name =='GSE122843':
    #     adata =  ut.specific_process(adata,dataname=data_name)
    # elif data_name =='GSE110894':
    #     adata =  ut.specific_process(adata,dataname=data_name)
    # elif data_name =='GSE112274':
    #     adata =  ut.specific_process(adata,dataname=data_name)
    # elif data_name =='GSE116237':
    #     adata =  ut.specific_process(adata,dataname=data_name)
    # elif data_name =='GSE108383':
    #     adata =  ut.specific_process(adata,dataname=data_name)
    # elif data_name =='GSE129730':
    #     adata =  ut.specific_process(adata,dataname=data_name)
    if data_name == 'GSE140440':
        adata = ut.specific_process(adata, dataname=data_name)
    elif data_name == 'GSE149383':
        adata = ut.specific_process(adata, dataname=data_name)
    else:
        adata = adata

    # Filter cells and genes
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    print(f"[Check] After filtering: cells = {adata.n_obs}, genes = {adata.n_vars}")

    adata = pp.cal_ncount_ngenes(adata)

    # Standardize the mitochondrial percentage column name
    if "pct_counts_mt" not in adata.obs.columns:
        if "percent.mt" in adata.obs.columns:
            adata.obs["pct_counts_mt"] = pd.to_numeric(
                adata.obs["percent.mt"],
                errors="coerce"
            )
        else:
            raise KeyError(
                "Neither 'pct_counts_mt' nor 'percent.mt' exists in metadata"
            )

    print(f"[Check] After calculating nCount/nGenes: cells = {adata.n_obs}, genes = {adata.n_vars}")

    print("Check cell count:", adata.n_obs)
    print("Check whether pct_counts_mt / n_genes_by_counts exist:",
          'pct_counts_mt' in adata.obs.columns,
          'n_genes_by_counts' in adata.obs.columns)
    print("Check data types:",
          adata.obs['pct_counts_mt'].dtype,
          adata.obs['n_genes_by_counts'].dtype)
    print("Check missing values:",
          adata.obs['pct_counts_mt'].isna().sum(),
          adata.obs['n_genes_by_counts'].isna().sum())
    print("Check minimum and maximum values:",
          adata.obs['pct_counts_mt'].min(), adata.obs['pct_counts_mt'].max(),
          adata.obs['n_genes_by_counts'].min(), adata.obs['n_genes_by_counts'].max())

    print("[DEBUG] Cell count before receipe_my:", adata.n_obs)
    sc.pp.filter_cells(adata, min_genes=200)
    print("[DEBUG] Cell count after filter_cells:", adata.n_obs)
    sc.pp.filter_genes(adata, min_cells=3)
    print("[DEBUG] Cell count after filter_genes:", adata.n_obs)

    # Preprocess data by filtering
    if data_name not in ['GSE112274', 'GSE140440']:
        adata = pp.receipe_my(adata, l_n_genes=min_n_genes, r_n_genes=max_n_genes, filter_mincells=args.min_c,
                              percent_mito=args.percent_mito,
                              filter_mingenes=args.min_g, normalize=True, log=True)
        print("[Main check] args.percent_mito =", args.percent_mito)
    else:
        adata = pp.receipe_my(adata, l_n_genes=min_n_genes, r_n_genes=max_n_genes, filter_mincells=args.min_c,
                              percent_mito=args.percent_mito, filter_mingenes=args.min_g,
                              normalize=True, log=True)
        print("[Main check] args.percent_mito =", args.percent_mito)

    print(f"[Check] After receipe_my preprocessing: cells = {adata.n_obs}, genes = {adata.n_vars}")

    # Select highly variable genes
    sc.pp.highly_variable_genes(adata, min_disp=g_disperson, max_disp=np.inf, max_mean=6)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]

    print(f"[Check] After highly variable gene selection: cells = {adata.n_obs}, genes = {adata.n_vars}")
    print(f"[Check] Number of highly variable genes = {adata.var.highly_variable.sum()}")


    # Preprocess data if spcific process is required
    data = adata.X
    print(f"[Check] Data matrix shape: {data.shape}")
    # PCA
    # Generate neighbor graph
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=10)
    # Generate cluster labels
    sc.tl.leiden(adata, resolution=leiden_res)
    sc.tl.umap(adata)
    adata.obs['leiden_origin'] = adata.obs['leiden']
    adata.obsm['X_umap_origin'] = adata.obsm['X_umap']
    data_c = adata.obs['leiden'].astype("long").to_list()

    ################################################# END SECTION OF SINGLE CELL DATA REPROCESSING ####################################################

    ################################################# START SECTION OF LOADING SC DATA TO THE TENSORS #################################################
    # Prepare to normailize and split target data
    mmscaler = preprocessing.MinMaxScaler()
    try:
        data = mmscaler.fit_transform(data)
    except:
        logging.warning("Only one class, no ROC")
        # Process sparse data
        data = data.todense()
        data = mmscaler.fit_transform(data)

    # Split data to train and valid set
    # Along with the leiden conditions for CVAE propose
    Xtarget_train, Xtarget_valid, Ctarget_train, Ctarget_valid = train_test_split(data, data_c, test_size=valid_size,
                                                                                  random_state=42)

    # Select the device of gpu
    if (args.device == "gpu"):
        device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
        torch.cuda.set_device(device)
        print('🚀 Using GPU ...')
    else:
        device = 'cpu'
        print('Using CPU ...')
    # Assuming that we are on a CUDA machine, this should print a CUDA device:
    logging.info(device)
    # Construct datasets and data loaders
    Xtarget_trainTensor = torch.FloatTensor(Xtarget_train).to(device)
    Xtarget_validTensor = torch.FloatTensor(Xtarget_valid).to(device)
    # print(Xtarget_validTensor.shape)
    # Use leiden label if CVAE is applied
    Ctarget_trainTensor = torch.LongTensor(Ctarget_train).to(device)
    Ctarget_validTensor = torch.LongTensor(Ctarget_valid).to(device)
    # print("C",Ctarget_validTensor )
    X_allTensor = torch.FloatTensor(data).to(device)
    C_allTensor = torch.LongTensor(data_c).to(device)

    train_dataset = TensorDataset(Xtarget_trainTensor, Ctarget_trainTensor)
    valid_dataset = TensorDataset(Xtarget_validTensor, Ctarget_validTensor)

    Xtarget_trainDataLoader = DataLoader(dataset=train_dataset, batch_size=batch_size, shuffle=True)
    Xtarget_validDataLoader = DataLoader(dataset=valid_dataset, batch_size=batch_size, shuffle=True)

    dataloaders_pretrain = {'train': Xtarget_trainDataLoader, 'val': Xtarget_validDataLoader}
    # print('START SECTION OF LOADING SC DATA TO THE TENSORS')
    ################################################# START SECTION OF LOADING SC DATA TO THE TENSORS #################################################

    ################################################# START SECTION OF LOADING BULK DATA  #################################################

    # Read source data
    data_r = pd.read_csv(source_data_path, index_col=0)
    label_r = pd.read_csv(label_path, index_col=0)
    if args.bulk == 'old':
        data_r = data_r[0:805]
        label_r = label_r[0:805]
    elif args.bulk == 'new':
        data_r = data_r[805:data_r.shape[0]]
        label_r = label_r[805:label_r.shape[0]]
    else:
        print("two databases combine")
    label_r = label_r.fillna(na)

    # Extract labels
    selected_idx = label_r.loc[:, select_drug] != na
    label = label_r.loc[selected_idx, select_drug]
    data_r = data_r.loc[selected_idx, :]
    label = label.values.reshape(-1, 1)

    # Encode labels
    le = preprocessing.LabelEncoder()
    label = le.fit_transform(label)
    dim_model_out = 2

    # Process source data
    mmscaler = preprocessing.MinMaxScaler()
    source_data = mmscaler.fit_transform(data_r)

    # Split source data
    Xsource_train_all, Xsource_test, Ysource_train_all, Ysource_test = train_test_split(source_data, label,
                                                                                        test_size=test_size,
                                                                                        random_state=42)
    Xsource_train, Xsource_valid, Ysource_train, Ysource_valid = train_test_split(Xsource_train_all, Ysource_train_all,
                                                                                  test_size=valid_size, random_state=42)

    # Transform source data
    # Construct datasets and data loaders
    Xsource_trainTensor = torch.FloatTensor(Xsource_train).to(device)
    Xsource_validTensor = torch.FloatTensor(Xsource_valid).to(device)

    Ysource_trainTensor = torch.LongTensor(Ysource_train).to(device)
    Ysource_validTensor = torch.LongTensor(Ysource_valid).to(device)

    sourcetrain_dataset = TensorDataset(Xsource_trainTensor, Ysource_trainTensor)
    sourcevalid_dataset = TensorDataset(Xsource_validTensor, Ysource_validTensor)

    Xsource_trainDataLoader = DataLoader(dataset=sourcetrain_dataset, batch_size=batch_size, shuffle=True)
    Xsource_validDataLoader = DataLoader(dataset=sourcevalid_dataset, batch_size=batch_size, shuffle=True)

    dataloaders_source = {'train': Xsource_trainDataLoader, 'val': Xsource_validDataLoader}
    # print('END SECTION OF LOADING BULK DATA')
    ################################################# END SECTION OF LOADING BULK DATA  #################################################

    ################################################# START SECTION OF MODEL CUNSTRUCTION  #################################################
    # Construct target encoder
    if reduce_model == "AE":
        encoder = AEBase(input_dim=data.shape[1], latent_dim=dim_au_out, h_dims=encoder_hdims, drop_out=args.dropout)
        loss_function_e = nn.MSELoss()
    elif reduce_model == "VAE":
        encoder = VAEBase(input_dim=data.shape[1], latent_dim=dim_au_out, h_dims=encoder_hdims, drop_out=args.dropout)
    if reduce_model == "DAE":
        encoder = AEBase(input_dim=data.shape[1], latent_dim=dim_au_out, h_dims=encoder_hdims, drop_out=args.dropout)
        loss_function_e = nn.MSELoss()

    logging.info("Target encoder structure is: ")
    logging.info(encoder)

    encoder.to(device)
    optimizer_e = optim.Adam(encoder.parameters(), lr=1e-2)
    loss_function_e = nn.MSELoss()
    exp_lr_scheduler_e = lr_scheduler.ReduceLROnPlateau(optimizer_e)

    # Binary classification
    dim_model_out = 2
    # Load the trained source encoder and predictor
    if reduce_model == "AE":
        source_model = PretrainedPredictor(input_dim=Xsource_train.shape[1], latent_dim=dim_au_out,
                                           h_dims=encoder_hdims,
                                           hidden_dims_predictor=predict_hdims, output_dim=dim_model_out,
                                           pretrained_weights=None, freezed=freeze, drop_out=args.dropout,
                                           drop_out_predictor=args.dropout)

        source_model.load_state_dict(torch.load(source_model_path))
        source_encoder = source_model
    if reduce_model == "DAE":
        source_model = PretrainedPredictor(input_dim=Xsource_train.shape[1], latent_dim=dim_au_out,
                                           h_dims=encoder_hdims,
                                           hidden_dims_predictor=predict_hdims, output_dim=dim_model_out,
                                           pretrained_weights=None, freezed=freeze, drop_out=args.dropout,
                                           drop_out_predictor=args.dropout)

        source_model.load_state_dict(torch.load(source_model_path))
        source_encoder = source_model
        # Load VAE model
    elif reduce_model in ["VAE"]:
        source_model = PretrainedVAEPredictor(input_dim=Xsource_train.shape[1], latent_dim=dim_au_out,
                                              h_dims=encoder_hdims,
                                              hidden_dims_predictor=predict_hdims, output_dim=dim_model_out,
                                              pretrained_weights=None, freezed=freeze, z_reparam=bool(args.VAErepram),
                                              drop_out=args.dropout, drop_out_predictor=args.dropout)
        source_model.load_state_dict(torch.load(source_model_path))
        source_encoder = source_model
    # logging.info("Load pretrained source model from: "+source_model_path)

    source_encoder.to(device)
    ################################################# END SECTION OF MODEL CUNSTRUCTION  #################################################

    ################################################# START SECTION OF SC MODEL PRETRAININIG  #################################################
    # Pretrain target encoder training
    # Pretain using autoencoder is pretrain is not False
    print("Pretrain training...")
    if (str(args.sc_encoder_path) != 'False'):
        # Pretrained target encoder if there are not stored files in the harddisk
        train_flag = True
        sc_encoder_path = str(sc_encoder_path)
        print("Pretrain==" + sc_encoder_path)

        # If pretrain is not False load from check point
        if (args.checkpoint != "False"):
            # if checkpoint is not False, load the pretrained model
            try:
                encoder.load_state_dict(torch.load(sc_encoder_path))
                logging.info("Load pretrained target encoder from " + sc_encoder_path)
                train_flag = False

            except:
                logging.info("Loading failed, procceed to re-train model")
                train_flag = True

        # If pretrain is not False and checkpoint is False, retrain the model
        if train_flag == True:

            if reduce_model == "AE":
                encoder, loss_report_en = t.train_AE_model(net=encoder, data_loaders=dataloaders_pretrain,
                                                           optimizer=optimizer_e, loss_function=loss_function_e,
                                                           load=False,
                                                           n_epochs=epochs, scheduler=exp_lr_scheduler_e,
                                                           save_path=sc_encoder_path)
            if reduce_model == "DAE":
                encoder, loss_report_en = t.train_DAE_model(net=encoder, data_loaders=dataloaders_pretrain,
                                                            optimizer=optimizer_e, loss_function=loss_function_e,
                                                            load=False,
                                                            n_epochs=epochs, scheduler=exp_lr_scheduler_e,
                                                            save_path=sc_encoder_path)

            elif reduce_model == "VAE":
                encoder, loss_report_en = t.train_VAE_model(net=encoder, data_loaders=dataloaders_pretrain,
                                                            optimizer=optimizer_e, load=False,
                                                            n_epochs=epochs, scheduler=exp_lr_scheduler_e,
                                                            save_path=sc_encoder_path)
            # print(loss_report_en)
            logging.info("Pretrained finished")

        # Before Transfer learning, we test the performance of using no transfer performance:
        # Use vae result to predict
        embeddings_pretrain = encoder.encode(X_allTensor)
        # print(embeddings_pretrain)
        pretrain_prob_prediction = source_model.predict(embeddings_pretrain).detach().cpu().numpy()
        adata.obs["sens_preds_pret"] = pretrain_prob_prediction[:, 1]
        adata.obs["sens_label_pret"] = pretrain_prob_prediction.argmax(axis=1)

        # Add embeddings to the adata object
        embeddings_pretrain = embeddings_pretrain.detach().cpu().numpy()
        adata.obsm["X_pre"] = embeddings_pretrain
    ################################################# END SECTION OF SC MODEL PRETRAININIG  #################################################

    ################################################# START SECTION OF TRANSFER LEARNING TRAINING #################################################
    # Using DaNN transfer learning
    # DaNN model
    # Set predictor loss
    loss_d = nn.CrossEntropyLoss()
    optimizer_d = optim.Adam(encoder.parameters(), lr=1e-2)
    exp_lr_scheduler_d = lr_scheduler.ReduceLROnPlateau(optimizer_d)

    # Set DaNN model
    # DaNN_model = DaNN(source_model=source_encoder,target_model=encoder)
    DaNN_model = DaNN(source_model=source_encoder, target_model=encoder, fix_source=bool(args.fix_source))
    DaNN_model.to(device)

    # Set distribution loss
    def loss(x, y, GAMMA=args.mmd_GAMMA):
        result = mmd.mmd_loss(x, y, GAMMA)
        return result

    loss_disrtibution = loss

    # Train DaNN model
    # logging.info("Trainig using" + mod + "model")
    target_model = TargetModel(source_model, encoder)
    # Switch to use regularized DaNN model or not
    print("DaNN Training...")
    if mod == 'ori':
        if args.checkpoint == 'True':
            DaNN_model, report_ = t.train_DaNN_model(DaNN_model,
                                                     dataloaders_source, dataloaders_pretrain,
                                                     # Should here be all optimizer d?
                                                     optimizer_d, loss_d,
                                                     epochs, exp_lr_scheduler_d,
                                                     dist_loss=loss_disrtibution,
                                                     load=target_model_path + "_DaNN.pkl",
                                                     weight=args.mmd_weight,
                                                     save_path=target_model_path + "_DaNN.pkl")
        else:
            DaNN_model, report_ = t.train_DaNN_model(DaNN_model,
                                                     dataloaders_source, dataloaders_pretrain,
                                                     # Should here be all optimizer d?
                                                     optimizer_d, loss_d,
                                                     epochs, exp_lr_scheduler_d,
                                                     dist_loss=loss_disrtibution,
                                                     load=False,
                                                     weight=args.mmd_weight,
                                                     save_path=target_model_path + "_DaNN.pkl")
    # Train DaNN model with new loss function
    if mod == 'new':
        # args.checkpoint = 'False'
        if args.checkpoint == 'True':
            DaNN_model, report_, _, _ = t.train_DaNN_model2(DaNN_model,
                                                            dataloaders_source, dataloaders_pretrain,
                                                            # Should here be all optimizer d?
                                                            optimizer_d, loss_d,
                                                            epochs, exp_lr_scheduler_d,
                                                            dist_loss=loss_disrtibution,
                                                            load=True,
                                                            weight=args.mmd_weight,
                                                            save_path=target_model_path + "_DaNN.pkl")
        else:
            DaNN_model, report_, _, _ = t.train_DaNN_model2(DaNN_model,
                                                            dataloaders_source, dataloaders_pretrain,
                                                            # Should here be all optimizer d?
                                                            optimizer_d, loss_d,
                                                            epochs, exp_lr_scheduler_d,
                                                            dist_loss=loss_disrtibution,
                                                            load=False,
                                                            weight=args.mmd_weight,
                                                            save_path=target_model_path + "_DaNN.pkl",
                                                            device=device)

    encoder = DaNN_model.target_model
    source_model = DaNN_model.source_model
    print("Transfer DaNN finished")
    ################################################# END SECTION OF TRANSER LEARNING TRAINING #################################################

    ################################################# START SECTION OF PREPROCESSING FEATURES #################################################
    # Extract feature embeddings
    # Extract prediction probabilities
    embedding_tensors = encoder.encode(X_allTensor)
    prediction_tensors = source_model.predictor(embedding_tensors)
    embeddings = embedding_tensors.detach().cpu().numpy()
    predictions = prediction_tensors.detach().cpu().numpy()
    print("predictions", predictions.shape)
    # Transform predict8ion probabilities to 0-1 labels

    adata.obs["sens_preds"] = predictions[:, 1]
    adata.obs["sens_label"] = predictions.argmax(axis=1)
    adata.obs["sens_label"] = adata.obs["sens_label"].astype('category')
    adata.obs["rest_preds"] = predictions[:, 0]

    ################################################# END SECTION OF ANALYSIS AND POST PROCESSING #################################################

    ################################################# START SECTION OF ANALYSIS FOR scRNA-Seq DATA #################################################
    # Save adata
    # adata.write("save/adata/"+data_name+para+".h5ad")
    ################################################# END SECTION OF ANALYSIS FOR scRNA-Seq DATA #################################################
    from sklearn.metrics import (average_precision_score,
                                 classification_report, mean_squared_error, r2_score, roc_auc_score, accuracy_score)
    report_df = {}
    Y_test_str = adata.obs['condition']
    Y_test = LabelEncoder().fit_transform(Y_test_str)

    sens_pb_results = adata.obs['sens_preds']
    lb_results = adata.obs['sens_label']

    # Y_test ture label
    ap_score = average_precision_score(Y_test, sens_pb_results)
    auroc_score = roc_auc_score(Y_test, sens_pb_results)
    print("auroc_score", auroc_score)
    print("ap_score", ap_score)
    acc = accuracy_score(Y_test, lb_results)
    print("Accuracy:", acc)
    report_dict = classification_report(Y_test, lb_results, output_dict=True)
    f1score = report_dict['weighted avg']['f1-score']
    precision = report_dict['weighted avg']['precision']
    recall = report_dict['weighted avg']['recall']
    print("precision", precision)
    print("recall", recall)
    print("f1score", f1score)
    report_df['f1_score'] = f1score
    baseline = np.mean(Y_test)
    norm_auprc = (ap_score - baseline) / (1 - baseline)
    print(f"Normalized AUPRC: {norm_auprc:.4f}")

    # file = 'save/bulk_f'+data_name+'_f1_score_ori.txt'
    # with open(file, 'a+') as f:
    #      f.write(para+'\t'+str(f1score)+'\n')
    print("sc model finished")
    sample_name = os.path.basename(args.sc_expression_file)
    sample_name = sample_name.replace("_raw_count_matrix.csv", "")

    result_filename = f"{sample_name}_{args.dimreduce}_bottle{args.bottleneck}_mod{args.mod}_{args.sampling}_prediction_results.csv"
    result_path = os.path.join(root_path, result_filename)
    adata.obs.to_csv(result_path)
    print(f"Prediction results saved to: {result_path}")
    # If print gene is true, then print gene
    # if (args.printgene=='T'):

    #     # Set up the TargetModel
    #     target_model = TargetModel(source_model,encoder)
    #     sc_X_allTensor=X_allTensor

    #     ytarget_allPred = target_model(sc_X_allTensor).detach().cpu().numpy()
    #     ytarget_allPred = ytarget_allPred.argmax(axis=1)
    #     # Caculate integrated gradient
    #     ig = IntegratedGradients(target_model)
    #     scattr, delta =  ig.attribute(sc_X_allTensor,target=1, return_convergence_delta=True,internal_batch_size=batch_size)
    #     scattr = scattr.detach().cpu().numpy()

    #     # Save integrated gradient
    #     igadata= sc.AnnData(scattr)
    #     igadata.var.index = adata.var.index
    #     igadata.obs.index = adata.obs.index
    #     sc_gra = "save/" + data_name +"sc_gradient.txt"
    #     sc_gen = "save/" + data_name +"sc_gene.csv"
    #     sc_lab = "save/" + data_name +"sc_label.csv"
    #     np.savetxt(sc_gra,scattr,delimiter = " ")
    #     DataFrame(adata.var.index).to_csv(sc_gen)
    #     DataFrame(adata.obs["sens_label"]).to_csv(sc_lab)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--sc_expression_file', type=str, required=True,
                        help='Path to the single cell expression data file (csv or other formats supported by pp.read_sc_file)')
    parser.add_argument('--meta_file', type=str, required=True,
                        help='Path to the metadata file')

    # data
    parser.add_argument('--raw_data_path', type=str, default='.', help='Raw data directory. Default: current directory')
    parser.add_argument('--root_path', type=str, default='.',
                        help='Root directory to store all outputs (logs, models, results). Default: current directory')
    parser.add_argument('--sc_data', type=str, default="GSE158457",
                        help='Accession id for testing data, only support pre-built data.')
    parser.add_argument('--meta', type=str, default=None, help="Path to metadata file (CSV format)")
    parser.add_argument('--drug', type=str, default='Dasatinib',
                        help='Name of the selected drug, should be a column name in the input file of --label')
    parser.add_argument('--checkpoint', type=str, default='False',
                        help='Load weight from checkpoint files, can be True,False, or a file path. Checkpoint files can be paraName1_para1_paraName2_para2... Default: True')
    parser.add_argument('--lr', type=float, default=0.01, help='Learning rate of model training. Default: 1e-2')
    parser.add_argument('--bottleneck', type=int, default=64,
                        help='Size of the bottleneck layer of the model. Default: 32')
    parser.add_argument('--dimreduce', type=str, default="DAE",
                        help='Encoder model type. Can be AE or VAE. Default: AE')
    parser.add_argument('--bulk_h_dims', type=str, default="512,256", help='Shape of the source encoder. Each number represent the number of neuron in a layer. \
                        Layers are seperated by a comma. Default: 512,256')
    parser.add_argument('--predictor_h_dims', type=str, default="256,128", help='Shape of the predictor. Each number represent the number of neuron in a layer. \
                        Layers are seperated by a comma. Default: 16,8')
    parser.add_argument('--dropout', type=float, default=0.1, help='Dropout of neural network. Default: 0.3')
    parser.add_argument('--sampling', type=str, default='upsampling', help='Samping method of training data for the bulk model traning. \
                        Can be no, upsampling, downsampling, or SMOTE. default: no')
    parser.add_argument('--device', type=str, default="gpu",
                        help='Device to train the model. Can be cpu or gpu. Deafult: cpu')

    parser.add_argument('--bulk_data', type=str, default='data/ALL_expression.csv',
                        help='Path of the bulk RNA-Seq expression profile')
    parser.add_argument('--label', type=str, default='data/ALL_label_binary_wf.csv',
                        help='Path of the processed bulk RNA-Seq drug screening annotation')
    parser.add_argument('--missing_value', type=int, default=1,
                        help='The value filled in the missing entry in the drug screening annotation, default: 1')
    parser.add_argument('--test_size', type=float, default=0.2,
                        help='Size of the test set for the bulk model traning, default: 0.2')
    parser.add_argument('--valid_size', type=float, default=0.2,
                        help='Size of the validation set for the bulk model traning, default: 0.2')
    parser.add_argument('--var_genes_disp', type=float, default=0, help='Dispersion of highly variable genes selection when pre-processing the data. \
                         If None, all genes will be selected .default: None')
    parser.add_argument('--min_n_genes', type=int, default=0,
                        help="Minimum number of genes for a cell that have UMI counts >1 for filtering propose, default: 0 ")
    parser.add_argument('--max_n_genes', type=int, default=20000,
                        help="Maximum number of genes for a cell that have UMI counts >1 for filtering propose, default: 20000 ")
    parser.add_argument('--min_g', type=int, default=200,
                        help="Minimum number of genes for a cell >1 for filtering propose, default: 200")
    parser.add_argument('--min_c', type=int, default=3,
                        help="Minimum number of cell that each gene express for filtering propose, default: 3")
    parser.add_argument('--percent_mito', type=int, default=100,
                        help="Percentage of expreesion level of moticondrial genes of a cell for filtering propose, default: 100")

    parser.add_argument('--cluster_res', type=float, default=0.2,
                        help="Resolution of Leiden clustering of scRNA-Seq data, default: 0.3")
    parser.add_argument('--mmd_weight', type=float, default=0.25,
                        help="Weight of the MMD loss of the transfer learning, default: 0.25")
    parser.add_argument('--mmd_GAMMA', type=int, default=1000,
                        help="Gamma parameter in the kernel of the MMD loss of the transfer learning, default: 1000")

    # trai
    parser.add_argument('--bulk_model_path', '-s', type=str, default='save/bulk_pre/',
                        help='Path of the trained predictor in the bulk level')
    parser.add_argument('--sc_model_path', '-p', type=str, default='save/sc_pre/',
                        help='Path (prefix) of the trained predictor in the single cell level')
    parser.add_argument('--sc_encoder_path', type=str, default='save/sc_encoder/',
                        help='Path of the pre-trained encoder in the single-cell level')
    parser.add_argument('--epochs', type=int, default=500, help='Number of epoches training. Default: 500')
    parser.add_argument('--batch_size', type=int, default=200, help='Number of batch size when training. Default: 200')
    parser.add_argument('--freeze_pretrain', type=int, default=0,
                        help='Fix the prarmeters in the pretrained model. 0: do not freeze, 1: freeze. Default: 0')
    parser.add_argument('--sc_h_dims', type=str, default="512,256", help='Shape of the encoder. Each number represent the number of neuron in a layer. \
                        Layers are seperated by a comma. Default: 512,256')
    parser.add_argument('--VAErepram', type=int, default=1)
    parser.add_argument('--batch_id', type=str, default="HN137", help="Batch id only for testing")
    parser.add_argument('--load_sc_model', type=int, default=0,
                        help='Load a trained model or not. 0: do not load, 1: load. Default: 0')
    parser.add_argument('--mod', type=str, default='new',
                        help='Embed the cell type label to regularized the training: new: add cell type info, ori: do not add cell type info. Default: new')
    parser.add_argument('--printgene', type=str, default='F', help='Print the cirtical gene list: T: print. Default: T')
    # miss
    parser.add_argument('--logging_file', '-l', type=str, default='save/logs/', help='Path of training log')
    parser.add_argument('--fix_source', type=int, default=0, help='Fix the bulk level model. Default: 0')
    parser.add_argument('--bulk', type=str, default='integrate',
                        help='Selection of the bulk database.integrate:both dataset. old: GDSC. new: CCLE. Default: integrate')
    #
    args, unknown = parser.parse_known_args()
    run_main(args)
