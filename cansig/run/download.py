"""Script for data downloading."""
import argparse

from typing import Optional


# Download choices
PIPELINE = "PIPELINE"
PREPROCESSING = "PREPROCESSING"


def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "what",
        help="The type of data to be downloaded. (Different tutorials use different datasets).",
        type=str,
        choices=[PIPELINE, PREPROCESSING],
        default=PIPELINE,
    )
    parser.add_argument(
        "--destination",
        help="Directory to which the data will be downloaded. If not specified, it will be created automatically.",
        type=str,
        default=None,
    )

    return parser


def download_hdf5(destination: Optional[str]) -> str:
    """This function downloads a small PIPELINE dataset for the pipeline
    demonstration.

    Args:
        destination: where the dataset should be downloaded


    Returns:
        the path to which the dataset was downloaded
        (if `destination` is not None, it should be `destination`).

    Note:
        The purpose of this function is to have a side effect.
    """
    # TODO(Pawel): Missing function.
    return "TODO"


def download_preprocessing(destination: Optional[str]) -> str:
    """This function downloads a dataset for the preprocessing tutorial.

    Args:
        destination: where the dataset should be downloaded


    Returns:
        the path to which the dataset was downloaded
        (if `destination` is not None, it should be `destination`).

    Note:
        The purpose of this function is to have a side effect.
    """
    # TODO(Pawel): Missing function.
    return "TODO"


def main() -> None:
    parser = create_parser()
    args = parser.parse_args()

    if args.what == PIPELINE:
        downloaded_path = download_hdf5(destination=args.destination)
    elif args.what == PREPROCESSING:
        downloaded_path = download_preprocessing(destination=args.destination)
    else:
        raise ValueError(f"Dataset {args.what} not recognized.")

    print(f"Download finished. The data is in {downloaded_path}.")


if __name__ == "__main__":
    main()
