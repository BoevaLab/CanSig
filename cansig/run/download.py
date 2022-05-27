"""Script for data downloading."""
import argparse
import pathlib
import shutil
import tempfile
from typing import Optional

import cansig.tutorial as tutorial

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
    destination = destination or "data/tutorial"
    destination = pathlib.Path(destination)
    destination.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmpdir:
        zip_path = pathlib.Path(tmpdir) / "data.zip"

        print("Downloading the ZIP version of the dataset...")

        tutorial.download_file_to_path(
            source_url="https://polybox.ethz.ch/index.php/s/NdjJLiNGcEwiswd/download",
            md5sum="7ea79c0e799b2b560c92cb2062f8d257",
            path=zip_path,
            verbose=True,
        )

        print("Downloaded the ZIP version. Unpacking...")
        shutil.unpack_archive(zip_path, extract_dir=destination)

    return str(destination)


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
