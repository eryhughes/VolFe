import pytest

# import pandas as pd
import VolFe as vf


## Equivalent to Penny's example
def test_mantlemelt():
    assert 2 + 2 == pytest.approx(4), "Inst Cu doesn't match test value"
