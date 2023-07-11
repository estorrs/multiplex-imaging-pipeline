from setuptools import setup, find_packages
from setuptools.command.install import install
from os import path
import subprocess

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md')) as f:
    long_description = f.read()

setup(
    # $ pip install multiplex-imaging-pipeline
    name='multiplex-imaging-pipeline',
    version='0.0.8',
    description='A Python library for multiplex imaging analysis',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/estorrs/multiplex-imaging-analysis',
    author='Erik Storrs',
    author_email='estorrs@wustl.edu',
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    keywords='multiplex imaging codex neighborhood analysis image segmentation visualization mibi codex phenocycler mihc hyperion',
    packages=find_packages(),
    python_requires='>=3.8',
    install_requires=[
        'pyyaml',
        'numpy',
        'pandas',
        'einops',
        'matplotlib',
        'scanpy',
        'seaborn',
        'tifffile',
        'ome-types',
        'scikit-image',
        'scikit-learn',
        'imagecodecs>=2022.7.27',
        'torch',
        'torchvision',
    ],
    extras_require={
        'segmentation': [
            'deepcell'
        ],
    },
    include_package_data=True,

    entry_points={
        'console_scripts': [
            'mip=multiplex_imaging_pipeline.multiplex_imaging_pipeline:main',
        ],
    },
)
