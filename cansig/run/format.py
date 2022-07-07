import argparse

import pandas as pd  # pytype: disable=import-error
import anndata as ad  # pytype: disable=import-error


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--expression-file", type=str, help="Path to the file containing expression data.")
    parser.add_argument("--observation-file", type=str, help="Name of the column with batch (or sample) index.")
    parser.add_argument("--filename", type=str, help="Name of the output file (with or without .h5ad suffix).")

    return parser


def output_path(raw: str) -> str:
    if raw.endswith(".h5ad"):
        return raw
    else:
        return f"{raw}.h5ad"


def main() -> None:
    # Read the CLI arguments
    parser = create_parser()
    args = parser.parse_args()

    expression_data = pd.read_csv(args.expression_file, index_col=0, sep=None)
    observation_data = pd.read_csv(args.observation_file, index_col=0, sep=None)

    adata = ad.AnnData(expression_data, obs=observation_data)

    adata.write_h5ad(output_path(args.filename), as_dense="X")
