#!/usr/bin/env python3

import setuptools

setuptools.setup(
    name="spaceranger-2-h5ad",
    version="0.0.1",
    author="Alma Andersson",
    author_email="almaan@kth.se",
    description="Converts spacreanger output to bundled h5ad file",
    url="https://github.com/almaan/space2h5ad",
    install_requires = [
        "pillow",
        "anndata",
        "h5py>=2.0.0",
        "pandas",
        "numpy",
        "scipy",
        ],
    packages = ['space2h5ad'],
    entry_points = {'console_scripts': ['space2h5ad = space2h5ad.__main__:main']},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.5',
)
