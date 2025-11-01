from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

long_description = (here / "README.md").read_text(encoding="utf-8")

setup(
    name="EIGENTOOLS",
    version="0.1.0",
    description="Python package designed for reading and writing PACKEDANCESTRYMAP files",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/floutt/EIGENTOOLS",
    author="Oluwatosin Olayinka",
    classifiers=[  # Optional
        "Development Status :: 4 - Beta",
        "Intended Audience :: Geneticists",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3 :: Only",
    ],
    keywords="EIGENSTRAT, PACKEDANCESTRYMAP, genetics, bioinformatics",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    python_requires=">=3.11",
    project_urls={
        "Bug Reports": "https://github.com/floutt/EIGENTOOLS/issues",
        "Source": "https://github.com/floutt/EIGENTOOLS",
    },
)
