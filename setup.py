# -*- coding: utf-8 -*-

# Extracting requirements from requirements.txt
import os
import io
# from Cython.Build import cythonize
from setuptools import setup, find_packages

with open('requirements.txt') as f:
    requirements = f.read().splitlines()


##########################
# File path read command
##########################
def read(*paths):
    """Build a file path from *paths* and return the contents."""
    with io.open(os.path.join(*paths), 'r', encoding='utf-8') as f:
        return f.read()


def setup_package():
    setup(

        name='MuPPI_dataset',

        version="0.0",
        python_requires='==3.8.5',
        packages=find_packages(),

        author="Dev Lis",
        author_email="contact.dev@lis-lab.fr",

        description="Code of MuPPI dataset generation",
        long_description=read('README.md'),

        include_package_data=True,
        url='https://gitlab.lis-lab.fr/dev/muppi-dataset-neurips/',
        install_requires=[requirements],

        classifiers=[
            "Programming Language :: Python",
            "Development Status :: 1 - Planning",
            "License :: OSI Approved",
            "Natural Language :: English",
            "Operating System :: Linux",
            "Programming Language :: Python :: 2/3",
            "Topic :: Machine Learning"],

        license="New BSD"

    )


if __name__ == "__main__":
    setup_package()
