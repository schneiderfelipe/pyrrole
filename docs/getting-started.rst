***************
Getting started
***************

In simple terms, the basic usage of pyrrole can be outlined in three steps:

1. Create a `data` object (this is actually a `pandas.DataFrame`).
2. Create a `ChemicalSystem` object.
3. Manipulate a `ChemicalSystem` object.

In order to understand each of them, let's walk through core API concepts as we tackle one everyday use case: the calculation of solvation free energy of acetic acid in water.

Solubility of acetic acid
=========================

Let's say that, after optimization and frequency calculations of `acetic acid <https://en.wikipedia.org/wiki/Acetic_acid>`_ were done (both in vacuuo and using an implicit solvation method [#level-of-theory]_), we wanted to calculate the `solvation energy <https://goldbook.iupac.org/html/S/ST07102.html>`_ of acetic acid in water.
This simple model perfectly exemplifies the usage of pyrrole, starting with the creation of a `data` object.

The `data` object
-----------------

The `data` object consists of a `pandas.DataFrame` whose records represent chemical species.
For our specific problem, we read logfiles (using the `read_cclib` function, which parses logfiles with the `cclib library <https://cclib.github.io/>`_) and store them in the required tabular form (using `create_data`):

>>> from pyrrole.atoms import read_cclib, create_data
>>> aa_vacuum = read_cclib("data/acetic_acid.out", "AcOH(g)")
>>> aa_water = read_cclib("data/acetic_acid@water.out", "AcOH(aq)")
>>> data = create_data([aa_vacuum, aa_water])

Each row of `data` above contains information found in a single logfile:

>>> columns = ["enthalpy", "entropy", "freeenergy"]
>>> data[columns]  # doctest: +NORMALIZE_WHITESPACE
            enthalpy  entropy  freeenergy
name
AcOH(g)  -228.533374 0.031135 -228.564509
AcOH(aq) -228.544332 0.030936 -228.575268

The energy values above are in `Hartree <https://en.wikipedia.org/wiki/Hartree>`_, which is the convention in the cclib project.
Learn more about `data` objects in `using-data-sets`.

The `ChemicalSystem` object
---------------------------

We are now in position to define our chemical system with `ChemicalSystem`.
Our model consists of a single equilibrium between gas phase and aqueous acetic acid:

>>> from pyrrole import ChemicalSystem
>>> system = ChemicalSystem("AcOH(g) <=> AcOH(aq)", data)
>>> system
ChemicalSystem(["AcOH(g) <=> AcOH(aq)"])

Usage of `ChemicalSystem`
~~~~~~~~~~~~~~~~~~~~~~~~~

`ChemicalSystem` objects can be manipulated in a variety of ways.
For instance, they can be converted to `pandas.DataFrame` (with the `ChemicalSystem.to_dataframe` method):

>>> reactions = system.to_dataframe()
>>> reactions[columns]  # doctest: +NORMALIZE_WHITESPACE
                      enthalpy   entropy  freeenergy
chemical_equation
AcOH(g) <=> AcOH(aq) -0.010958 -0.000198   -0.010759

Again, energy values are given in Hartree.
Conversion factors can be used for handling other units (with the help of the `scipy.constants` module):

>>> from scipy.constants import kilo, N_A, physical_constants
>>> hartree, _, _ = physical_constants["Hartree energy"]
>>> factor = hartree * N_A / kilo  # Hartree to kJ/mol
>>> factor
2625.4996382852164

The calculated factor can be used to convert a whole table if so desired:

>>> reactions[columns] * factor  # doctest: +NORMALIZE_WHITESPACE
                      enthalpy   entropy  freeenergy
chemical_equation
AcOH(g) <=> AcOH(aq) -28.76991 -0.521109  -28.248775

(By the way, the reported experimental value for the solvation free energy of acetic acid in water is -28.0 kJ/mol [#experimental-freeenergy-acetic-acid]_.)

Now we're ready to start `using-data-sets`.

.. [#level-of-theory] Calculations were done at `PBEh-3c`_/`SMD`_ (water) using the `ORCA`_ electronic structure package (version 4.0.1.2). Logfiles can be found in the `project's repository <https://github.com/dudektria/pyrrole>`_.

.. _`PBEh-3c`: https://doi.org/10.1063/1.4927476
.. _`SMD`: https://doi.org/10.1021/jp810292n
.. _`ORCA`: https://orcaforum.cec.mpg.de/

.. [#experimental-freeenergy-acetic-acid] *J. Phys. Chem. B*, **2009**, *113* (18), pp 6378-6396 DOI: `10.1021/jp810292n <https://doi.org/10.1021/jp810292n>`_ (supporting information).
