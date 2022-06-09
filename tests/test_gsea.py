from typing import Tuple

import pytest  # pytype: disable=import-error

import cansig.gsea as cgsea


@pytest.mark.parametrize("input_output", [("SOME_PATHWAY", "SOME\nPATHWAY"), ("SOME PATHWAY", "SOME\nPATHWAY")])
def test_default_formatter(input_output: Tuple[str, str]) -> None:
    inp, out = input_output

    formatter = cgsea.DefaultFormatter()
    out_ = formatter.format(inp)
    assert out_ == out
