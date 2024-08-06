[![](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)

Welcome to VolFe! an open-source framework for calculating melt-vapour equilibria including silicate melt, carbon, hydrogen, sulfur and noble gases. 

This work is currently in preparation for publication. 

For more information, see the ReadTheDocs page (more content will be added soon)!

https://volfe.readthedocs.io/en/latest/

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