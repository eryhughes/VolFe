from setuptools import setup, find_packages
from os import path

this_directory = path.abspath(path.dirname(__file__))

with open(
    path.join(this_directory, "src", "VolFe", "_version.py"), encoding="utf-8"
) as f:
    exec(f.read())

# with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
#    long_description = f.read()


setup(
    name="VolFe",
    version=__version__,
    author="Ery Hughes",
    author_email="ery.c.hughes@gmail.com",
    description="VolFe",
    # long_description=long_description,
    # long_description_content_type="text/markdown",
    url="https://github.com/eryhughes/VolFe",
    package_dir={"": "src"},  # Optional
    packages=find_packages(where="src"),  # Required
    # package_data={
    #    # Include all pickle files
    #    "": ["*.pkl"],
    # },
    install_requires=[
        "pandas",
        "numpy<2",
        "datetime",
        "gmpy2",
        "scipy",
        "densityx",
        "PySulfSat",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    extras_require={
        "dev": ["pytest", "sympy"],
    },
)
