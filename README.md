# VolFe

Welcome to VolFe! An open-source framework for calculating melt-vapour equilibria including silicate melt, carbon, hydrogen, sulfur, and noble gases.

[![](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![PyPI](https://badgen.net/pypi/v/VolFe)](https://pypi.org/project/VolFe/)
[![Build Status](https://github.com/eryhughes/VolFe/actions/workflows/main.yml/badge.svg?branch=main)](https://github.com/eryhughes/VolFe/actions/workflows/main.yml)
[![Documentation Status](https://readthedocs.org/projects/VolFe/badge/?version=latest)](https://VolFe.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/646746387.svg)](https://doi.org/10.5281/zenodo.15756761)

Read more about the python package in the Volcanica article for the chemistry and speciation calculations (and please cite if you use the package!):

Hughes EC, Liggins P, Wieser P, and Stolper EM (2025). VolFe: an open-source Python package for calculating melt-vapor equilibria including silicate melt, carbon, hydrogen, sulfur, and noble gases. Volcanica 8(2):457-481 https://doi.org/10.30909/vol/imvc1781

The Supplementary Material of the Lithos article explains the isotope calculations (and please additionally cite this paper if you use the isotope calculations):

Saper LM, Bromiley G, Cao R, Brounce M, Hughes EC, and Woelki D (2025). The primary magmatic δ34S of the Troodos Ophiolite and evidence for early and late sulfide saturation. Lithos 518-519:108331 https://doi.org/10.1016/j.lithos.2025.108331


For more information and worked examples, see the ReadTheDocs page:
https://volfe.readthedocs.io/en/latest/

VolFe can be installed using pip from PyPI or from GitHub (see notes below about installing an editable version). Please install VolFe in its own environment because it will rewrite all your package versions.

## Development

If you wish to edit VolFe on your own computer, you can install an editable version using

```
pip install -e ".[dev]"
```
from inside a virtual environment (use either venv or anaconda). This will import VolFe
in a format that allows you to run any edits you have made, and all it's requirements,
alongside useful packages for developing VolFe (pytest, sympy).

Check VolFe runs on your machine, and that any edits you make haven't broken existing code by running pytest:
```
python -m pytest tests
```
or you can use the testing frameworks that come with your IDE (e.g. [VSCode](https://code.visualstudio.com/docs/python/testing), [PyCharm](https://www.jetbrains.com/help/pycharm/testing-your-first-python-application.html)).
