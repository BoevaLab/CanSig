import pydantic  # pytype: disable=import-error
import pytest  # pytype: disable=import-error

import cansig.filesys as fs


@pytest.mark.parametrize("mode", ("file", "dir"))
def test_get_file__ok(tmp_path, mode: str) -> None:
    file = tmp_path / "a.txt"
    file.touch()

    if mode == "file":
        assert fs.get_file(file_or_dir=file, ext=".txt") == file
    elif mode == "dir":
        assert fs.get_file(file_or_dir=tmp_path, ext=".txt") == file


def test_get_file__doesnt_exist(tmp_path) -> None:
    with pytest.raises(FileNotFoundError):
        assert fs.get_file(tmp_path, "txt")

    with pytest.raises(FileNotFoundError):
        assert fs.get_file(tmp_path / "doesnt_exist_dir", "")


@pytest.mark.parametrize("n_files", (2, 5))
def test_get_file__many(tmp_path, n_files: int) -> None:
    for i in range(n_files):
        (tmp_path / f"{i}.txt").touch()
    with pytest.raises(FileExistsError):
        fs.get_file(tmp_path, ".txt")


def test_saving_and_loading(tmp_path) -> None:
    class Grocery(pydantic.BaseModel):
        apples: int
        milk: float

    path = tmp_path / "target"
    # Save settings
    settings = Grocery(apples=3, milk=2.0)
    fs.save_settings(settings=settings, path=path)
    # Read them and compare
    settings_ = fs.read_settings(factory=Grocery, path=path)
    assert settings_ == settings
