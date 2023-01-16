"""The high-level utilities for the integration step."""
import argparse
import pathlib
import logging
from typing import cast, Protocol, Union

import anndata  # pytype: disable=import-error
import pandas as pd  # pytype: disable=import-error

import cansig.filesys as fs  # pytype: disable=import-error
import cansig.logger as clogger  # pytype: disable=import-error
import cansig.models.api as models  # pytype: disable=import-error

LOGGER = logging.getLogger(__name__)
DEFAULT_OUTPUT_BASE_PATH = pathlib.Path("./outputs/batch-integration")


class Arguments(Protocol):
    """Protocol used to define CLI arguments."""

    @property
    def model(self) -> str:
        """Integration model to be used."""
        raise NotImplementedError

    @property
    def path(self) -> pathlib.Path:
        raise NotImplementedError

    @property
    def batch(self) -> str:
        """Batch column name."""
        raise NotImplementedError

    @property
    def data(self) -> pathlib.Path:
        """The path to the AnnData object with malignant cells."""
        raise NotImplementedError

    @property
    def latent(self) -> int:
        """The dimension of the latent space."""
        raise NotImplementedError

    @property
    def output(self) -> pathlib.Path:
        """Output directory."""
        raise NotImplementedError

    @property
    def max_epochs(self) -> int:
        """Maximum number of epochs to be run."""
        raise NotImplementedError

    @property
    def log(self) -> pathlib.Path:
        """Location of the log file."""
        raise NotImplementedError

    @property
    def n_top_genes(self) -> int:
        """Number of most highly variable genes to be used."""
        raise NotImplementedError


def parse_args() -> Arguments:
    """Creates the CLI parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument("data", type=pathlib.Path, help="The path to the anndata object.")
    parser.add_argument("batch", type=str, help="Name of the batch column.")
    parser.add_argument("--latent", type=int, help="The dimensionality of the latent space.", default=10)
    parser.add_argument("--max-epochs", type=int, help="The maximal number of training epochs.", default=400)
    parser.add_argument("--model", type=str, default="scvi", choices=["scvi", "cansig"])
    default_output = DEFAULT_OUTPUT_BASE_PATH / fs.get_directory_name()
    parser.add_argument("--output", type=pathlib.Path, help="Output directory.", default=default_output)
    parser.add_argument(
        "--log", type=pathlib.Path, help="Where the log file should be saved.", default=pathlib.Path("integration.log")
    )
    parser.add_argument(
        "--n-top-genes", type=int, default=2000, help="The number of most highly variable genes to use."
    )

    args = parser.parse_args()
    return cast(Arguments, args)


def integrate_adata(
    data: anndata.AnnData,
    config: Union[models.SCVIConfig, models.CanSigConfig],
) -> pd.DataFrame:
    if isinstance(config, models.SCVIConfig):
        model = models.SCVI(config=config, data=data)
        representations = model.get_latent_codes()
    elif isinstance(config, models.CanSigConfig):
        model = models.CanSigWrapper(config=config, data=data)
        representations = model.get_latent_codes()
    else:
        raise NotImplementedError(f"Something went wrong, unknown config {config}.")
    return representations


def integrate(
    data_path: pathlib.Path,
    config: models.SCVIConfig,
    output: pathlib.Path,
) -> bool:
    # Save settings
    output_dir = fs.IntegrationDir(output, create=True)
    fs.save_settings(settings=config, path=output_dir.integration_settings)

    # Train the model and get the representations
    data = anndata.read_h5ad(data_path)
    representations = integrate_adata(data=data, config=config)

    # Save the representations
    fs.save_latent_representations(representations=representations, path=output_dir.latent_representations)

    return output_dir.valid()


def main(args: Arguments) -> None:
    clogger.configure_logging(args.log)

    LOGGER.info(f"Initializing integration model config for {args.model}...")
    # pytype: disable=attribute-error
    if args.model == "scvi":
        config = models.SCVIConfig(
            batch=args.batch,
            n_latent=args.latent,
            preprocessing=models.module_scvi.PreprocessingConfig(n_top_genes=args.n_top_genes),
        )
    elif args.model == "cansig":
        config = models.CanSigConfig(
            batch=args.batch,
            n_latent=args.latent,
            preprocessing=models.module_cansig.PreprocessingConfig(n_top_genes=args.n_top_genes),
        )
    else:
        raise NotImplementedError(f"Model {args.model} is not implemented.")

    # pytype: enable=attribute-error
    # Set the number of training epochs.
    # It's a bit hacky, as we could do that at the initialization stage.
    config.train.max_epochs = args.max_epochs
    integrate(data_path=args.data, config=config, output=args.output)

    LOGGER.info("Run finished.")


if __name__ == "__main__":
    main(parse_args())
