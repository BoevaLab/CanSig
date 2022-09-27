from typing import Tuple

import pydantic  # pytype: disable=import-error
import pytest  # pytype: disable=import-error

import cansig.gsea as cgsea  # pytype: disable=import-error


@pytest.mark.parametrize("input_output", [("SOME_PATHWAY", "SOME\nPATHWAY"), ("SOME PATHWAY", "SOME\nPATHWAY")])
def test_default_formatter(input_output: Tuple[str, str]) -> None:
    inp, out = input_output

    formatter = cgsea.DefaultFormatter()
    out_ = formatter.format(inp)
    assert out_ == out


def test_invalid_gmt(tmp_path):
    path = tmp_path / "wrong.gmt"
    with open(path, "w") as f:
        f.write("This is a wrong GMT file!")

    with pytest.raises(pydantic.ValidationError):
        cgsea.GeneExpressionConfig(gene_sets=path)


def test_gmt_doesnt_exist():
    with pytest.raises(pydantic.ValidationError):
        cgsea.GeneExpressionConfig(gene_sets="this-file-doesnt-exist.gmt")
