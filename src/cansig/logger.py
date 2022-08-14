"""Setting up logger. Based on Hydra configuration."""
import logging.config
import pathlib
from typing import Union

_Filename = Union[str, pathlib.Path]


def default_config(filename: _Filename) -> dict:
    """Generates a dictionary config for `logging.config.dictConfig`.

    The logs will be streamed to the standard output and saved to `filename`."""
    filename = str(filename)

    return {
        "version": 1,
        "formatters": {"simple": {"format": "[%(asctime)s][%(name)s][%(levelname)s] - %(message)s"}},
        "handlers": {
            "console": {
                "class": "logging.StreamHandler",
                "formatter": "simple",
                "stream": "ext://sys.stdout",
            },
            "file": {
                "class": "logging.FileHandler",
                "formatter": "simple",
                "filename": filename,
            },
        },
        "root": {
            "level": "INFO",
            "handlers": ["console", "file"],
        },
        "disable_existing_loggers": False,
    }


def configure_logging(filename: _Filename) -> None:
    """Configures logging. Logs will be streamed to the standard output and saved to `filename`."""
    conf = default_config(filename)
    logging.config.dictConfig(conf)
