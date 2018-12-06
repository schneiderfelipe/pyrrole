pyrrole
=======

.. |pypi-badge| image:: https://badge.fury.io/py/pyrrole.svg
   :target: https://badge.fury.io/py/pyrrole
   :alt: Python Package Index (PyPI)

.. |build-badge| image:: https://travis-ci.org/dudektria/pyrrole.svg?branch=master
   :target: https://travis-ci.org/dudektria/pyrrole
   :alt: Travis CI Status

.. |docs-badge| image:: https://img.shields.io/badge/docs-pyrrole-blue.svg
   :target: https://pyrrole.readthedocs.io/en/latest/?badge=latest
   :alt: Latest Documentation

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
You could also want to read tabular data :ref:`from a file <Reading local files>` (or even :ref:`from the web <Reading the web>`) using `pandas <https://pandas.pydata.org/>`_.

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

Interested? `Have another example <getting-started>`.

.. [#level-of-theory] Optimizations and frequency calculations of both ammonia and the planar transition state were performed at `PBEh-3c`_ using the `ORCA`_ electronic structure package (version 4.0.1.2). Logfiles can be found in the `project's repository <https://github.com/dudektria/pyrrole/tree/master/data>`_.

.. _`PBEh-3c`: https://doi.org/10.1063/1.4927476
.. _`ORCA`: https://orcaforum.cec.mpg.de/

.. [#experimental-freeenergy-ammonia-inversion] *Chem. Phys. Lett.*, **2003**, *370* (3), pp 360-365 DOI: `10.1016/S0009-2614(03)00107-6 <https://doi.org/10.1016/S0009-2614(03)00107-6>`_.

Installation
------------

You can get the library directly from `PyPI <https://pypi.org/project/pyrrole/>`_:

.. code-block:: console

   $ pip install pyrrole

.. image:: data:image/svg+xml;base64,PD94bWwgdmVyc2lvbj0iMS4wIj8+CjxzdmcgdmVyc2lvbj0iMS4xIiBpZD0idG9wc3ZnIgp4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHhtbG5zOnhsaW5rPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hsaW5rIgp4bWxuczpjbWw9Imh0dHA6Ly93d3cueG1sLWNtbC5vcmcvc2NoZW1hIiB4PSIwIiB5PSIwIiB3aWR0aD0iMjAwcHgiIGhlaWdodD0iMjAwcHgiIHZpZXdCb3g9IjAgMCAxMDAgMTAwIj4KPHRpdGxlPnB5cnJvbGUgLSBPcGVuIEJhYmVsIERlcGljdGlvbjwvdGl0bGU+CjxyZWN0IHg9IjAiIHk9IjAiIHdpZHRoPSIxMDAiIGhlaWdodD0iMTAwIiBmaWxsPSJ3aGl0ZSIvPgo8ZyB0cmFuc2Zvcm09InRyYW5zbGF0ZSgwLDApIj4KPHN2ZyB3aWR0aD0iMTAwIiBoZWlnaHQ9IjEwMCIgeD0iMCIgeT0iMCIgdmlld0JveD0iMCAwIDE0NC43MjEgMTQxLjU1NCIKZm9udC1mYW1pbHk9InNhbnMtc2VyaWYiIHN0cm9rZT0icmdiKDAsMCwwKSIgc3Ryb2tlLXdpZHRoPSIyIiAgc3Ryb2tlLWxpbmVjYXA9InJvdW5kIj4KPGxpbmUgeDE9IjkyLjQiIHkxPSIxMDEuNiIgeDI9IjEwNC43IiB5Mj0iNjMuNSIgc3Ryb2tlPSJyZ2IoMCwwLDApIiAgc3Ryb2tlLXdpZHRoPSIyLjAiLz4KPGxpbmUgeDE9Ijg3LjQiIHkxPSI5My42IiB4Mj0iOTYuMCIgeTI9IjY3LjAiIHN0cm9rZT0icmdiKDAsMCwwKSIgIHN0cm9rZS13aWR0aD0iMi4wIi8+CjxsaW5lIHgxPSIxMDQuNyIgeTE9IjYzLjUiIHgyPSI4Mi45IiB5Mj0iNDcuNiIgc3Ryb2tlPSJyZ2IoMCwwLDApIiAgc3Ryb2tlLXdpZHRoPSIyLjAiLz4KPGxpbmUgeDE9IjYxLjgiIHkxPSI0Ny42IiB4Mj0iNDAuMCIgeTI9IjYzLjUiIHN0cm9rZT0icmdiKDAsMCwwKSIgIHN0cm9rZS13aWR0aD0iMi4wIi8+CjxsaW5lIHgxPSI0MC4wIiB5MT0iNjMuNSIgeDI9IjUyLjQiIHkyPSIxMDEuNiIgc3Ryb2tlPSJyZ2IoMCwwLDApIiAgc3Ryb2tlLXdpZHRoPSIyLjAiLz4KPGxpbmUgeDE9IjQ4LjciIHkxPSI2Ny4wIiB4Mj0iNTcuNCIgeTI9IjkzLjYiIHN0cm9rZT0icmdiKDAsMCwwKSIgIHN0cm9rZS13aWR0aD0iMi4wIi8+CjxsaW5lIHgxPSI1Mi40IiB5MT0iMTAxLjYiIHgyPSI5Mi40IiB5Mj0iMTAxLjYiIHN0cm9rZT0icmdiKDAsMCwwKSIgIHN0cm9rZS13aWR0aD0iMi4wIi8+Cjx0ZXh0IHg9IjY2LjM2MDY4MCIgeT0iNDguMDAwMDAwIiBmaWxsPSJyZ2IoMTIsMTIsMjU1KSIgIHN0cm9rZT0icmdiKDEyLDEyLDI1NSkiIHN0cm9rZS13aWR0aD0iMSIgZm9udC1zaXplPSIxNiIgPk48L3RleHQ+Cjx0ZXh0IHg9IjY2LjM2MDY4MCIgeT0iMzIuMDAwMDAwIiBmaWxsPSJyZ2IoMTIsMTIsMjU1KSIgIHN0cm9rZT0icmdiKDEyLDEyLDI1NSkiIHN0cm9rZS13aWR0aD0iMSIgZm9udC1zaXplPSIxNiIgPkg8L3RleHQ+Cjwvc3ZnPgo8Y21sOm1vbGVjdWxlIGlkPSJweXJyb2xlIj4KIDxjbWw6YXRvbUFycmF5PgogIDxjbWw6YXRvbSBpZD0iYTEiIGVsZW1lbnRUeXBlPSJDIiB4Mj0iMC44MDkwMTciIHkyPSItMC41ODc3ODUiLz4KICA8Y21sOmF0b20gaWQ9ImEyIiBlbGVtZW50VHlwZT0iQyIgeDI9IjAuNTAwMDAwIiB5Mj0iLTEuNTM4ODQyIi8+CiAgPGNtbDphdG9tIGlkPSJhMyIgZWxlbWVudFR5cGU9IkMiIHgyPSItMC41MDAwMDAiIHkyPSItMS41Mzg4NDIiLz4KICA8Y21sOmF0b20gaWQ9ImE0IiBlbGVtZW50VHlwZT0iQyIgeDI9Ii0wLjgwOTAxNyIgeTI9Ii0wLjU4Nzc4NSIvPgogIDxjbWw6YXRvbSBpZD0iYTUiIGVsZW1lbnRUeXBlPSJOIiB4Mj0iMC4wMDAwMDAiIHkyPSItMC4wMDAwMDAiLz4KIDwvY21sOmF0b21BcnJheT4KIDxjbWw6Ym9uZEFycmF5PgogIDxjbWw6Ym9uZCBhdG9tUmVmczI9ImEyIGEzIiBvcmRlcj0iMSIvPgogIDxjbWw6Ym9uZCBhdG9tUmVmczI9ImEyIGExIiBvcmRlcj0iMiIvPgogIDxjbWw6Ym9uZCBhdG9tUmVmczI9ImEzIGE0IiBvcmRlcj0iMiIvPgogIDxjbWw6Ym9uZCBhdG9tUmVmczI9ImExIGE1IiBvcmRlcj0iMSIvPgogIDxjbWw6Ym9uZCBhdG9tUmVmczI9ImE0IGE1IiBvcmRlcj0iMSIvPgogPC9jbWw6Ym9uZEFycmF5Pgo8L2NtbDptb2xlY3VsZT4KPC9nPgo8L3N2Zz4K
  :alt: pyrrole
  :align: center
