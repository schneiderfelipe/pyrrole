*********************
Systems and equations
*********************

The `ChemicalSystem` object revisited
=====================================

In `getting-started`, we saw the basics of chemical systems.

Drawing.

Internally, a `ChemicalSystem` object consists of individual `ChemicalEquation` objects, which can be manipulated on their own.

The `ChemicalEquation` object
=============================

Single chemical equations in pyrrole are handled by `ChemicalEquation` objects.
A special mini-language is used to define chemical equations in a way that
makes it easy to simply copy and paste from the web.
For instance, the following metal displacement was obtained from a `Wikipedia entry`_:

.. _Wikipedia entry: https://en.wikipedia.org/wiki/Redox#Metal_displacement

>>> from pyrrole import ChemicalEquation
>>> half_zinc = ChemicalEquation('Zn(s) -> Zn+2(aq) + 2 e-')
>>> half_copper = ChemicalEquation('Cu+2(aq) + 2 e- <- Cu(s)')

`ChemicalEquation` objects can be manipulated just like vectors, i.e., summed and multiplied by scalar values:

>>> half_zinc - half_copper
ChemicalEquation('Cu+2(aq) + Zn(s) -> Cu(s) + Zn+2(aq)')

Stoichiometry coefficients can be obtained individually:

>>> half_zinc.coefficient['e-']
2.0

There's no need to use chemical formulae for chemical species.
Any mix of printable characters can be used:

>>> ChemicalEquation('cis-A <=> trans-A')
ChemicalEquation('cis-A <=> trans-A')
