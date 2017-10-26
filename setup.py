#!/usr/bin/env python
from setuptools import setup, find_packages
import json
import glob

scripts = glob.glob('scripts/*.py')

if __name__ == '__main__':
    setup(
        packages=find_packages(),
        name="asetk",
        author="Leopold Talirz",
        author_email="leopold.talirz@gmail.com",
        description="Toolkit for working with atomic and electronic structure data, built on top of the Atomistic Simulation Environment (ASE) for atomic structures.",
        url="https://github.com/ltalirz/asetk",
        license="MIT",
        classifiers= [
            "Programming Language :: Python"
        ],
        version="0.3.0",
        install_requires=[
            "numpy >= 1.9",
            "scipy >= 0.14",
            "matplotlib >= 1.4",
            "ase >= 3.8.4",
        ],
        scripts=scripts,
    )
