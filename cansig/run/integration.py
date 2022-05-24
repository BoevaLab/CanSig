import argparse
import pathlib
from typing import cast, Protocol

import anndata  # pytype: disable=import-error

import cansig.filesys as fs
import cansig.models.api as models

DEFAULT_OUTPUT_BASE_PATH = pathlib.Path("./outputs/batch-integration")


class Arguments(Protocol):
    @property
    def n_latent(self) -> int:
        raise NotImplementedError

    @property
    def path(self) -> pathlib.Path:
        raise NotImplementedError

    @property
    def batch(self) -> str:
        raise NotImplementedError

    @property
    def data(self) -> pathlib.Path:
        raise NotImplementedError

    @property
    def latent(self) -> int:
        raise NotImplementedError

    @property
    def output(self) -> pathlib.Path:
        raise NotImplementedError

    @property
    def max_epochs(self) -> int:
        raise NotImplementedError


def parse_args() -> Arguments:
    parser = argparse.ArgumentParser()
    parser.add_argument("data", type=pathlib.Path, help="The path to the anndata object.")
    parser.add_argument("batch", type=str, help="Name of the batch column.")
    parser.add_argument("--latent", type=int, help="The dimensionality of the latent space.", default=10)
    parser.add_argument("--max-epochs", type=int, help="The maximal number of training epochs.", default=400)

    default_output = DEFAULT_OUTPUT_BASE_PATH / fs.get_directory_name()
    parser.add_argument("--output", type=pathlib.Path, help="Output directory.", default=default_output)

    args = parser.parse_args()
    return cast(Arguments, args)


def integrate(data_path: pathlib.Path, config: models.SCVIConfig, output: pathlib.Path,) -> bool:
    # Save settings
    output_dir = fs.IntegrationDir(output, create=True)
    fs.save_settings(settings=config, path=output_dir.integration_settings)

    # Train the model and get the representations
    data = anndata.read_h5ad(data_path)
    model = models.SCVI(config=config, data=data)
    representations = model.get_latent_codes()

    # Save the representations
    fs.save_latent_representations(representations=representations, path=output_dir.latent_representations)

    return output_dir.valid()


def main(args: Arguments) -> None:
    config = models.SCVIConfig(batch=args.batch, n_latent=args.latent,)
    # Set the number of training epochs.
    # It's a bit hacky, as we could do that at the initialization stage.
    config.train.max_epochs = args.max_epochs

    integrate(data_path=args.data, config=config, output=args.output)


if __name__ == "__main__":
    main(parse_args())
