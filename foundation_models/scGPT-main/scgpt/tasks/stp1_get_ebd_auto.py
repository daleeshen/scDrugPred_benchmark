import argparse
import json
import os
from pathlib import Path
from typing import Union

import anndata
import numpy as np
import pandas as pd
import pickle
import torch

from anndata import AnnData
from torch.utils.data import DataLoader, SequentialSampler
from tqdm import tqdm

from data_collator_EBD import DataCollator
from model import TransformerModel
from tokenizer import GeneVocab
from utils import load_pretrained

PathLike = Union[str, os.PathLike]


def get_batch_cell_embeddings(adata, cell_embedding_mode: str = "cls", model=None, vocab=None, max_length=1200,
                              batch_size=64, model_configs=None, gene_ids=None, use_batch_labels=False) -> np.ndarray:
    count_matrix = adata.X
    count_matrix = count_matrix if isinstance(count_matrix, np.ndarray) else count_matrix.A

    # gene vocabulary ids
    if gene_ids is None:
        gene_ids = np.array(adata.var["id_in_vocab"])
        assert np.all(gene_ids >= 0)

    class Dataset(torch.utils.data.Dataset):
        def __init__(self, count_matrix, gene_ids, batch_ids=None):
            self.count_matrix = count_matrix
            self.gene_ids = gene_ids
            self.batch_ids = batch_ids

        def __len__(self):
            return len(self.count_matrix)

        def __getitem__(self, idx):
            row = self.count_matrix[idx]
            nonzero_idx = np.nonzero(row)[0]
            values = row[nonzero_idx]
            genes = self.gene_ids[nonzero_idx]

            # append <cls> token at the beginning
            genes = np.insert(genes, 0, vocab["<cls>"])
            values = np.insert(values, 0, model_configs["pad_value"])
            genes = torch.from_numpy(genes).long()
            values = torch.from_numpy(values).float()

            output = {"id": idx, "genes": genes, "expressions": values}
            if self.batch_ids is not None:
                output["batch_labels"] = self.batch_ids[idx]
            return output

    if cell_embedding_mode == "cls":
        dataset = Dataset(count_matrix, gene_ids)
        collator = DataCollator(do_padding=True, pad_token_id=vocab[model_configs["pad_token"]],
                                pad_value=model_configs["pad_value"], do_mlm=False, do_binning=True,
                                max_length=max_length, sampling=True, keep_first_n_tokens=1)
        data_loader = DataLoader(dataset, batch_size=batch_size, sampler=SequentialSampler(dataset),
                                 collate_fn=collator, drop_last=False,
                                 num_workers=min(len(os.sched_getaffinity(0)), batch_size), pin_memory=True,
                                 shuffle=False)
        device = next(model.parameters()).device
        cell_embeddings = np.zeros((len(dataset), model_configs["embsize"]), dtype=np.float32)

        with torch.no_grad(), torch.cuda.amp.autocast(enabled=True):
            count = 0
            for data_dict in tqdm(data_loader, desc="Embedding cells"):
                input_gene_ids = data_dict["gene"].to(device)
                src_key_padding_mask = input_gene_ids.eq(vocab[model_configs["pad_token"]])
                embeddings = model._encode(input_gene_ids, data_dict["expr"].to(device),
                                           src_key_padding_mask=src_key_padding_mask,
                                           batch_labels=data_dict["batch_labels"].to(
                                               device) if use_batch_labels else None)

                # get CLS embedding
                embeddings = embeddings[:, 0, :].cpu().numpy()
                cell_embeddings[count: count + len(embeddings)] = embeddings
                count += len(embeddings)

        # L2 normalization
        cell_embeddings = cell_embeddings / np.linalg.norm(cell_embeddings, axis=1, keepdims=True)
    else:
        raise ValueError(f"Unknown cell embedding mode: {cell_embedding_mode}")

    return cell_embeddings


def embed_data(adata_or_file: Union[AnnData, PathLike], model_dir: PathLike, output_file: PathLike,
               gene_col: str = "index", max_length=1200, batch_size=64, use_fast_transformer=True):
    csv_file_path = Path(adata_or_file)
    print(f"Reading input: {csv_file_path}")

    # expression matrix
    df = pd.read_csv(
        csv_file_path,
        sep="\t"
    )

    cell_barcode = df["Cell_barcode"].astype(str).tolist()
    condition = df["Condition"].tolist()
    label = [1 if x == "sensitive" else 0 for x in condition]

    df = df.drop(columns=["Cell_barcode", "Condition"])

    # metadata and labels
    metadata = pd.read_csv(csv_file_path, sep="\t", usecols=["Cell_barcode", "Condition"])
    cell_barcode = metadata["Cell_barcode"].astype(str).tolist()
    condition = metadata["Condition"].tolist()
    label = [1 if element == "sensitive" else 0 for element in condition]

    print(f"Cells: {df.shape[0]}")
    print(f"Genes before vocab matching: {df.shape[1]}")
    print(f"Sensitive: {sum(label)}")
    print(f"Resistant: {len(label) - sum(label)}")

    # AnnData
    adata = anndata.AnnData(df)
    adata.X = adata.X.astype(np.float32)

    # gene names
    adata.var["index"] = adata.var.index
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    # ----------------------
    # Load pretrained model
    # ----------------------
    model_dir = Path(model_dir)
    vocab_file = model_dir / "vocab.json"
    model_config_file = model_dir / "args.json"
    model_file = model_dir / "best_model.pt"

    # Check model files
    for required_file in [vocab_file, model_config_file, model_file]:
        if not required_file.exists():
            raise FileNotFoundError(f"Model file not found: {required_file}")

    pad_token = "<pad>"
    special_tokens = [pad_token, "<cls>", "<eoc>"]

    # ----------------------
    # Vocabulary
    # ----------------------
    vocab = GeneVocab.from_file(vocab_file)
    for s in special_tokens:
        if s not in vocab:
            vocab.append_token(s)

    adata.var["id_in_vocab"] = [vocab[gene] if gene in vocab else -1 for gene in adata.var[gene_col]]

    # keep genes in vocabulary
    adata = adata[:, adata.var["id_in_vocab"] >= 0]
    print(f"Genes after vocab matching: {adata.shape[1]}")
    if adata.shape[1] == 0:
        raise ValueError("No genes matched the scGPT vocabulary.")

    # model config
    with open(model_config_file, "r") as f:
        model_configs = json.load(f)

    vocab.set_default_index(vocab["<pad>"])
    genes = adata.var[gene_col].tolist()
    gene_ids = np.array(vocab(genes), dtype=int)

    # ----------------------
    # Build model
    # ----------------------
    model = TransformerModel(ntoken=len(vocab), d_model=model_configs["embsize"], nhead=model_configs["nheads"],
                             d_hid=model_configs["d_hid"], nlayers=model_configs["nlayers"],
                             nlayers_cls=model_configs["n_layers_cls"], n_cls=1, vocab=vocab,
                             dropout=model_configs["dropout"], pad_token=model_configs["pad_token"],
                             pad_value=model_configs["pad_value"], do_mvc=True, do_dab=False, use_batch_labels=False,
                             domain_spec_batchnorm=False, explicit_zero_prob=False,
                             use_fast_transformer=use_fast_transformer, fast_transformer_backend="flash",
                             pre_norm=False)
    load_pretrained(model, torch.load(model_file, map_location=device), verbose=False)
    model.to(device)
    model.eval()

    # ----------------------
    # Generate embeddings
    # ----------------------
    cell_embeddings = get_batch_cell_embeddings(adata, cell_embedding_mode="cls", model=model, vocab=vocab,
                                                max_length=max_length, batch_size=batch_size,
                                                model_configs=model_configs, gene_ids=gene_ids, use_batch_labels=False)
    print("Embedding shape:", cell_embeddings.shape)

    # ----------------------
    # Save
    # ----------------------
    data_dict = {"embeds": cell_embeddings, "labels": label,
                 "cell_barcode": cell_barcode, "condition": condition}

    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "wb") as p:
        pickle.dump(data_dict, p)

    print(f"Saved embedding: {output_file}")


def main():
    parser = argparse.ArgumentParser(description="Generate scGPT cell embeddings")
    parser.add_argument("--input_file", type=str, required=True, help="Input txt file")
    parser.add_argument("--input_root", type=str, required=True, help="Root directory of input data, e.g. scdata/")
    parser.add_argument("--output_root", type=str, required=True,
                        help="Root directory for embeddings, e.g. scgpt/samples/")
    parser.add_argument("--model_dir", type=str, required=True,
                        help="Directory containing vocab.json, args.json and best_model.pt")
    parser.add_argument("--batch_size", type=int, default=64)
    parser.add_argument("--max_length", type=int, default=1200)
    parser.add_argument("--gene_col", type=str, default="index")
    args = parser.parse_args()

    input_file = Path(args.input_file).resolve()
    input_root = Path(args.input_root).resolve()
    output_root = Path(args.output_root).resolve()
    model_dir = Path(args.model_dir).resolve()

    # ----------------------
    # Check input
    # ----------------------
    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")
    if input_file.suffix != ".txt":
        raise ValueError(f"Input file must be txt: {input_file}")

    # ----------------------
    # Get relative path
    # ----------------------
    try:
        relative_path = input_file.relative_to(input_root)
    except ValueError:
        raise ValueError(f"\nInput file:\n{input_file}\n\nis not inside input root:\n{input_root}")

    # Example: imbalance/1_30/1/GSE108397.txt -> imbalance/1_30/1/GSE108397_samples.pkl
    output_dir = output_root / relative_path.parent
    output_file = output_dir / (relative_path.stem + "_samples.pkl")

    print("\n" + "=" * 70)
    print(f"Input file : {input_file}")
    print(f"Relative   : {relative_path}")
    print(f"Output file: {output_file}")
    print(f"Model dir  : {model_dir}")
    print("=" * 70 + "\n")

    embed_data(adata_or_file=input_file, model_dir=model_dir, output_file=output_file, use_fast_transformer=True,
               max_length=args.max_length, batch_size=args.batch_size, gene_col=args.gene_col)


if __name__ == "__main__":
    main()