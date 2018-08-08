********************
Using `data` objects
********************

Any `pandas.DataFrame` indexed by names of chemical species is a valid `data` object in pyrrole [#standard-gibbs-free-energy-of-formation]_:

>>> import pandas as pd
>>> data = pd.DataFrame(
...     [{'name': 'CO3-2(aq)', 'freeenergy': -527.8},
...      {'name': 'HCO3-(aq)', 'freeenergy': -586.85},
...      {'name': 'H2CO3(aq)', 'freeenergy': -623.1},
...      {'name': 'OH-(aq)', 'freeenergy': -157.2},
...      {'name': 'H2O(l)', 'freeenergy': -237.14}])
>>> data = data.set_index('name')
>>> data  # doctest: +NORMALIZE_WHITESPACE
           freeenergy
name
CO3-2(aq)     -527.80
HCO3-(aq)     -586.85
H2CO3(aq)     -623.10
OH-(aq)       -157.20
H2O(l)        -237.14

The `pandas library <https://pandas.pydata.org/>`_, a dependency of pyrrole, can be used to create `data` objects.
Below are examples of creating `data` objects from different sources.

Reading local files
===================

Pandas can read data sets in various formats, such as
`comma-separated values (CSV) <https://en.wikipedia.org/wiki/Comma-separated_values>`_,
`Google BigQuery <https://en.wikipedia.org/wiki/BigQuery>`_,
`Hierarchical Data Format (HDF) <https://en.wikipedia.org/wiki/Hierarchical_Data_Format>`_,
`JavaScript Object Notation (JSON) <http://www.json.org/>`_,
`Microsoft Excel <https://en.wikipedia.org/wiki/Microsoft_Excel>`_,
and many `other supported format types <https://pandas.pydata.org/pandas-docs/stable/io.html>`_:

>>> data = pd.read_hdf("data/data.h5")
>>> data[['jobfilename', 'freeenergy', 'enthalpy']]
                  jobfilename  freeenergy    enthalpy
0            data/acetate.out -228.000450 -227.969431
1      data/acetate@water.out -228.120113 -228.089465
2        data/acetic_acid.out -228.564509 -228.533374
3  data/acetic_acid@water.out -228.575268 -228.544332

Pyrrole requires indices to represent names of chemical species, which is, like above, not always the case.
Setting meaningful indices can be accomplished by feeding a custom function to `data.apply`:

>>> def update(series):
...     """Compute a new column 'name' and add it to row."""
...     series['name'] = (series['jobfilename']
...                       .replace('data/', '')
...                       .replace('.out', ''))
...     series['name'] = (series['name']
...                       .replace('acetate', 'AcO-')
...                       .replace('acetic_acid', 'AcOH'))
...     series['name'] = series['name'].replace('@water', '(aq)')
...     if '(aq)' not in series['name']:
...         series['name'] += "(g)"
...     return series

The function above should be applied to the `data` object, which can then be reindexed:

>>> data = data.apply(update, axis='columns').set_index('name')
>>> data[['jobfilename', 'freeenergy', 'enthalpy']]  # doctest: +NORMALIZE_WHITESPACE
                         jobfilename  freeenergy    enthalpy
name
AcO-(g)             data/acetate.out -228.000450 -227.969431
AcO-(aq)      data/acetate@water.out -228.120113 -228.089465
AcOH(g)         data/acetic_acid.out -228.564509 -228.533374
AcOH(aq)  data/acetic_acid@water.out -228.575268 -228.544332

The `data` object is now ready to be used:

>>> from pyrrole import ChemicalSystem
>>> system = ChemicalSystem(['AcO-(g) <=> AcO-(aq)',
...                          'AcOH(g) <=> AcOH(aq)'],
...                         data['freeenergy'])
>>> system.to_dataframe()  # doctest: +NORMALIZE_WHITESPACE
                      freeenergy
chemical_equation
AcO-(g) <=> AcO-(aq)   -0.119663
AcOH(g) <=> AcOH(aq)   -0.010759

In `getting-started`, we showed how to use `create_data` to produce a `data` object by reading output files from computational chemistry programs.
Reading lots of logfiles is slow, which is why storing the data in a file translates to faster retrievals later.
This can be accomplished with `ccframe <http://cclib.github.io/how_to_parse.html#ccframe>`_, a command-line tool that is part of `cclib <http://cclib.github.io/>`_ (a dependency of pyrrole).
In fact, the file ``data.h5`` used in the example above was produced using ccframe:

.. code-block:: console

   $ ccframe -O data/data.h5 data/acetate*out data/acetic_acid*out

Learn more about ccframe in both its help page (``$ ccframe -h``) and `documentation <http://cclib.github.io/how_to_parse.html#ccframe>`_.

Reading the web
===============

There's a lot of freely available data on the internet.
For instance, `NIST <https://www.nist.gov/>`_ offers `enthalpies of formation at 0K <https://cccbdb.nist.gov/hf0k.asp>`_ (in kJ/mol).
Luckily, pandas supports `reading HTML tables <https://pandas.pydata.org/pandas-docs/stable/io.html#html>`_ directly:

>>> url = "https://cccbdb.nist.gov/hf0k.asp"
>>> data = pd.read_html(url, header=0)[3]  # fourth table in page
>>> data = data.set_index("Species")
>>> data = data[["Name", "Hfg 0K", "DOI"]]
>>> data.head()  # doctest: +NORMALIZE_WHITESPACE
                         Name  Hfg 0K                       DOI
Species
D              Deuterium atom   219.8                       NaN
H               Hydrogen atom   216.0  10.1002/bbpc.19900940121
H+       Hydrogen atom cation  1528.1                       NaN
D2         Deuterium diatomic     0.0                       NaN
H2          Hydrogen diatomic     0.0  10.1002/bbpc.19900940121

This data allows us to calculate the `bond-dissociation enthalpy <https://en.wikipedia.org/wiki/Bond-dissociation_energy>`_ of the hydrogen molecule at 0K, for instance:

>>> from pyrrole import ChemicalEquation
>>> equation = ChemicalEquation("H2 -> 2 H", data)
>>> equation.to_series()
Hfg 0K    432.0
Name: H2 -> 2 H, dtype: float64

That's 432 kJ/mol, or 103.3 kcal/mol.

It's time to take a deeper look at `systems-and-equations`.

.. [#standard-gibbs-free-energy-of-formation] Obtained from `standard Gibbs free energy of formation <https://en.wikipedia.org/wiki/Standard_Gibbs_free_energy_of_formation>`_.
