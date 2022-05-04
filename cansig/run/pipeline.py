"""The main pipeline.

Takes as input the data and multi-run specification, and then processes the data according
to all models specified.
In the end, produces summary.
"""
import argparse
import logging
import pathlib

logger = logging.getLogger(__name__)


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("data", type=pathlib.Path, help="HDF5 file containing the dataset.")
    parser.add_argument(
        "--model-runs",
        type=int,
        default=1,
        help="Number of data integration models to be trained per the specified latent dimension. "
        "Each model will be trained with a different random seed.",
    )
    parser.add_argument("--cluster-runs", type=int, default=1, help="Number of random seeds used for the clustering.")
    parser.add_argument(
        "--dimensions", nargs="+", default=[4, 6], help="List with the number of latent dimensions to be used."
    )
    parser.add_argument(
        "--clusters", nargs="+", default=[3, 5, 10], help="List with the number of clusters to be used."
    )
    parser.add_argument(
        "--output", type=pathlib.Path, default=pathlib.Path("output"), help="Output directory."
    )  # TODO(Pawel): Parametrize this with datetime.
    parser.add_argument(
        "--save-intermediate",
        type=bool,
        default=True,
        help="Whether the intermediate results should be saved."
        "By default, the results are saved. Turning this off is discouraged, "
        "unless the system memory is very limited.",
    )
    return parser


def main() -> None:
    parser = create_parser()
    args = parser.parse_args()

    print(args)


if __name__ == "__main__":
    main()
