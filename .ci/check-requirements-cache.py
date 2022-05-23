"""Checks whether the requirements specified in
`setup.cfg` and `requirements-dev.txt` sum up to `requirements.txt`.

Rationale:
    - `setup.cfg`: dependencies needed to run the package
    - `requirements-dev.txt`: requirements needed to run the CI tools as well
    - `requirements.txt`: file used by GitHub for caching
"""
import configparser
import sys
from typing import List, Set


def parse_config_lines(raw: str) -> List[str]:
    """Returns the lines in a config section."""
    return [x for x in raw.split("\n") if x]


def read_install_requires(config: configparser.ConfigParser) -> List[str]:
    """Reads the `install_requires` section of `setup.cfg`."""
    raw = config["options"]["install_requires"]
    return parse_config_lines(raw)


def read_additional(config: configparser.ConfigParser) -> List[str]:
    """Reads the test requirements of `setup.cfg`."""
    raw = config["options.extras_require"]["test"]
    return parse_config_lines(raw)


def read_requirements(filepath) -> Set[str]:
    """Reads a standard requirements.txt file."""
    with open(filepath) as fh:
        raw = fh.readlines()
    raw = [x.strip() for x in raw]
    raw = [x for x in raw if x and x[0] != "#"]
    return set(raw)


def read_setup_dependencies(path) -> Set[str]:
    """Reads all the dependencies in `setup.cfg`."""
    config = configparser.ConfigParser()
    config.read(path)

    return set(read_install_requires(config) + read_additional(config))


def main():
    setup_deps = read_setup_dependencies("setup.cfg")
    requirements_dev = read_requirements("requirements-dev.txt")
    requirements_txt = read_requirements("requirements.txt")
    requirements_total = setup_deps.union(requirements_dev)

    if requirements_total != requirements_txt:
        diff = requirements_total.symmetric_difference(requirements_txt)
        print(f"Requirements differ: {diff}.")
        sys.exit(1)


if __name__ == "__main__":
    main()
