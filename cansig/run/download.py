"""Script for data downloading."""
import argparse

from typing import Optional


# Download choices
HDF5 = "HDF5"
PREPROCESSING = "PREPROCESSING"


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "what",
        help="The type of data to be downloaded. (Different tutorials use different datasets).",
        type=str,
        choices=[HDF5, PREPROCESSING],
        default=HDF5,
    )
    parser.add_argument(
        "--destination",
        help="Directory to which the data will be downloaded. If not specified, it will be created automatically.",
        type=str,
        default=None,
    )

    return parser


def download_hdf5(destination: Optional[str]) -> str:
    pass


def download_preprocessing(destination: Optional[str]) -> str:
    pass


def main() -> None:
    parser = create_parser()
    args = parser.parse_args()

    if args.what == HDF5:
        downloaded_path = download_hdf5(destination=args.destination)
    elif args.what == PREPROCESSING:
        downloaded_path = download_preprocessing(destination=args.destination)
    else:
        raise ValueError(f"Dataset {args.what} not recognized.")

    print(f"Download finished. The data is in {downloaded_path}.")


if __name__ == "__main__":
    main()
