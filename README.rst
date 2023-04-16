**WARNING**: this project has been retired.
Please migrate to
`geem-lab/overreact <https://github.com/geem-lab/overreact>`_.

pyrrole
=======

.. |pypi-badge| image:: https://badge.fury.io/py/pyrrole.svg
   :target: https://badge.fury.io/py/pyrrole
   :alt: Python Package Index (PyPI)

.. |build-badge| image:: https://travis-ci.org/dudektria/pyrrole.svg?branch=master
   :target: https://travis-ci.org/dudektria/pyrrole
   :alt: Travis CI Status

.. |docs-badge| image:: https://readthedocs.org/projects/pyrrole/badge/?version=latest
   :target: https://pyrrole.readthedocs.io/en/latest/?badge=latest
   :alt: Latest Documentation Status

|docs-badge| |build-badge| |pypi-badge|

A Python package for solving chemical problems with computational modeling.

Usage example
-------------

As a usage example, let's calculate the energy barrier involved in `nitrogen inversion in ammonia <https://en.wikipedia.org/wiki/Nitrogen_inversion>`_.

.. figure:: https://upload.wikimedia.org/wikipedia/commons/2/2d/Nitrogen-inversion-3D-balls.png
   :alt: Nitrogen inversion in ammonia
   :align: center

   When ammonia turns "inside out", it passes through a planar transition state (`image in public domain <https://commons.wikimedia.org/wiki/File:Nitrogen-inversion-3D-balls.png>`_).

We do this in three simple steps (only eight lines of code):

1. **Get the data**

We first obtain the raw data, which will later be fed to our chemical model.
Below we read computational chemistry logfiles of both ground and transition states [#level-of-theory]_.

>>> from pyrrole.atoms import read_cclib, create_data
>>> gs = read_cclib("data/ammonia/ammonia.out", name="NH3(g)")
>>> ts = read_cclib("data/ammonia/invers.out", name="NH3(g)#")
>>> data = create_data(gs, ts)

Pyrrole uses `cclib <https://cclib.github.io/>`_ for reading logfiles, which is `compatible with all major computational chemistry packages <https://cclib.github.io/#summary>`_.
You could also want to read tabular data `from a file <https://pyrrole.readthedocs.io/en/latest/using-data-sets.html#reading-local-files>`_ (or even `from the web <https://pyrrole.readthedocs.io/en/latest/using-data-sets.html#reading-the-web>`_) using `pandas <https://pandas.pydata.org/>`_.

2. **Specify the model**

We now describe our model.
This is accomplished through chemical equations:

>>> from pyrrole import ChemicalEquation
>>> equation = ChemicalEquation("NH3(g) -> NH3(g)#", data)

While model above consists of a single `ChemicalEquation`, you could create complex models with multiple chemical equations with `ChemicalSystem` objects.
You might want to store your complex models in text files too.

3. **Obtain the results**

Simply let pyrrole calculate the energy barrier:

>>> results = equation.to_series()
>>> results["freeenergy"] * 2625.4996382852164  # Hartree to kJ/mol
19.30952589472923

(As a side note, the reference value is 21.162 kJ/mol [#experimental-freeenergy-ammonia-inversion]_.)

Interested? `Have another example <https://pyrrole.readthedocs.io/en/latest/getting-started.html>`_.

.. [#level-of-theory] Optimizations and frequency calculations of both ammonia and the planar transition state were performed at `PBEh-3c`_ using the `ORCA`_ electronic structure package (version 4.0.1.2). Logfiles can be found in the `project's repository <https://github.com/dudektria/pyrrole/tree/master/data>`_.

.. _`PBEh-3c`: https://doi.org/10.1063/1.4927476
.. _`ORCA`: https://orcaforum.cec.mpg.de/

.. [#experimental-freeenergy-ammonia-inversion] *Chem. Phys. Lett.*, **2003**, *370* (3), pp 360-365 DOI: `10.1016/S0009-2614(03)00107-6 <https://doi.org/10.1016/S0009-2614(03)00107-6>`_.

Installation
------------

You can get the library directly from `PyPI <https://pypi.org/project/pyrrole/>`_:

.. code-block:: console

   $ pip install pyrrole
