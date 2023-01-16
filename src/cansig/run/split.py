"""Script used to split a data set by patients into two data sets."""
from typing import Tuple

import argparse

import anndata as ad  # pytype: disable=import-error
import numpy as np
import scanpy as sc  # pytype: disable=import-error


def create_parser() -> argparse.ArgumentParser:
    """Creates the CLI parser."""
    parser = argparse.ArgumentParser(description="Split a data set by patients into two data sets.")

    parser.add_argument("COLUMN", type=str, help="Column storing patient IDs, to be split by.")
    parser.add_argument("INPUT", help="Path to the .h5ad object with the data.")
    parser.add_argument("OUTPUT1", help="Output path used to save the first split.")
    parser.add_argument("OUTPUT2", help="Output path used to save the second split.")
    parser.add_argument(
        "--random-seed", type=int, default=152, help="Random seed used to split the data set into two groups."
    )

    return parser


def split_by_column(
    data: ad.AnnData,
    column: str,
    seed: int,
) -> Tuple[ad.AnnData, ad.AnnData]:
    """Splits a data set into two basing on the ``column`` values.

    Args:
        data: data set to be splitted
        column: column name storing categorical values
        seed: random seed for reproducibility

    Returns:
        anndata, a slice (not a copy) with one random split
        anndata, a slice (not a copy) with the other random split
    """

    if column not in data.obs.columns:
        raise ValueError(f"Column {column} not in columns: {data.obs.columns}.")

    all_values = list(data.obs[column].unique())

    if len(all_values) < 2:
        raise ValueError(f"Too few column values to split into two groups: {all_values}.")

    rng = np.random.default_rng(seed)
    rng.shuffle(all_values)

    split_index = len(all_values) // 2

    split1 = all_values[:split_index]
    split2 = all_values[split_index:]

    print(f"Splitting {len(data)} cells by {column} into two groups: {split1} and {split2}...")

    index1 = data.obs[column].isin(split1)
    index2 = data.obs[column].isin(split2)

    data1 = data[index1]
    data2 = data[index2]

    print(f"The first data set has {len(data1)} cells and the second has {len(data2)} cells.")

    return data1, data2


def main() -> None:
    """The main function.
    Parses the CLI arguments, generates the splits and saves them to H5AD files."""
    parser = create_parser()
    args = parser.parse_args()
    outputs = [args.OUTPUT1, args.OUTPUT2]

    print(f"Reading input file {args.INPUT}...")

    data = sc.read_h5ad(args.INPUT)

    print("Calculating data splits...")
    splits = split_by_column(data=data, column=args.COLUMN, seed=args.random_seed)

    print("Saving the data...")
    for split, output_file in zip(splits, outputs):
        split.write_h5ad(output_file)

    print("Finished!")


if __name__ == "__main__":
    main()
