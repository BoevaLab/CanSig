import pathlib
from typing import List

import cansig.filesys as fs


class MultirunDirectory(fs.StructuredDir):
    def valid(self) -> bool:
        # TODO(Pawel): Consider making this more elaborate.
        return True

    @property
    def integration_directories(self) -> pathlib.Path:
        return self.path / "integration"

    @property
    def postprocessing_directories(self) -> pathlib.Path:
        return self.path / "postprocessing"

    @property
    def metasig_directories(self) -> pathlib.Path:
        return self.path / "metasignatures"

    @property
    def analysis_directory(self) -> pathlib.Path:
        return self.path / "final_analysis"


def get_valid_dirs(multirun_dir: MultirunDirectory) -> List[fs.PostprocessingDir]:
    valid_dirs = []
    invalid_dirs = []
    for path in multirun_dir.postprocessing_directories.glob("*"):
        wrapped_path = fs.PostprocessingDir(path)
        if wrapped_path.valid():
            valid_dirs.append(wrapped_path)
        else:
            invalid_dirs.append(path)

    if invalid_dirs:
        print(f"The following directories seem to be corrupted: {invalid_dirs}.")

    return valid_dirs
