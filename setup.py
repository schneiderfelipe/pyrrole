#!/usr/bin/python3

"""Package setup script for pyrrole."""

import os.path
import setuptools
from distutils.core import setup

name = 'pyrrole'
this_directory = os.path.abspath(os.path.dirname(__file__))

version_file = open(os.path.join(this_directory, "VERSION"))
version = version_file.read().strip()

url = 'https://github.com/dudektria/pyrrole'
download_url = \
   '{:s}/archive/{:s}.tar.gz'.format(url, version)

with open(os.path.join(this_directory, 'README.rst')) as f:
    long_description = f.read()

setup(name=name,
      version=version,
      url=url,
      download_url=download_url,
      author='Felipe S. S. Schneider',
      author_email='schneider.felipe@posgrad.ufsc.br',
      license='MIT',
      description=('A Python package for solving problems in chemistry with '
                   'computational modeling.'),
      long_description=long_description,
      classifiers=["Programming Language :: Python :: 3",
                   "Programming Language :: Python",
                   "License :: OSI Approved :: MIT License",
                   "Operating System :: OS Independent",
                   "Development Status :: 3 - Alpha",
                   "Environment :: Console",
                   "Intended Audience :: Science/Research",
                   "Intended Audience :: Education",
                   "Intended Audience :: Developers",
                   "Topic :: Scientific/Engineering :: Chemistry",
                   "Topic :: Education",
                   "Topic :: Software Development :: Libraries :: Python Modules"],  # noqa
      keywords=['science',
                'research',
                'chemistry'],
      packages=setuptools.find_packages(exclude=['*test*']),
      install_requires=["cclib>=1.5.3",
                        "matplotlib>=2.1.1",
                        "networkx>=2.1",
                        "openbabel>=2.4.1",
                        "pandas>=0.23.4",
                        "pyparsing>=2.2.0",
                        "tables"],
      setup_requires=['nose>=1.0'],
      test_suite='nose.collector',
      include_package_data=True,
      command_options={
        'build_sphinx': {
            'project': ('setup.py', name),
            'version': ('setup.py', version),
            'release': ('setup.py', version),
            'build_dir': ('setup.py', 'docs/_build'),
            'builder': ('setup.py', 'html doctest'),
            }
      },
      )
