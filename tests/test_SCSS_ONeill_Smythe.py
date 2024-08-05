# Equivalent to Penny's example

# run with: python -m pytest tests
# And/or set up using VSCode if that's what you use

import pytest


def test_mantlemelt():
    assert 2 + 2 == pytest.approx(4), "Inst Cu doesn't match test value"
