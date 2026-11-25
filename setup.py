from setuptools import setup, find_packages
from pathlib import Path

# see https://packaging.python.org/guides/single-sourcing-package-version/
version_dict = {}
with open(Path(__file__).parents[0] / "chempotpy/_version.py") as fp:
    exec(fp.read(), version_dict)
version = version_dict["__version__"]
del version_dict

setup(
    name="chempotpy",
    version=version,
    author=(
        "Yinan Shu, Zoltan Varga, Dayou Zhang, Donald G. Truhlar, "
        "Ruben Meana-Paneda"
    ),
    description=(
        "Self-contained Python and modern Fortran library of potential "
        "energy surfaces"
    ),
    long_description=(Path(__file__).parent / "README.md").read_text(),
    long_description_content_type="text/markdown",
    url="https://github.com/rmeanapa/chempotpy",
    project_urls={
        "Documentation": "https://github.com/rmeanapa/chempotpy#readme",
        "Source": "https://github.com/rmeanapa/chempotpy",
        "Submit a PES": (
            "https://github.com/rmeanapa/chempotpy/issues/new"
            "?template=new-pes.yml"
        ),
        "Issues": "https://github.com/rmeanapa/chempotpy/issues",
    },
    license="Apache-2.0",
    python_requires=">=3.8",
    packages=find_packages(),
    package_data={
        "chempotpy": ["*/*.f90", "core/*.json"],
    },
    install_requires=[
        'numpy',
    ],
    zip_safe=False,
)
