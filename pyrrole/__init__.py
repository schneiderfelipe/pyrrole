#!/usr/bin/python3
# -*- encoding: utf-8 -*-

"""
Basic default fuctionality imported from submodules.

Both `ChemicalEquation` and `ChemicalSystem` are provided:

>>> from pyrrole import ChemicalEquation, ChemicalSystem
>>> ChemicalEquation('NaCl(s) <=> Na+(aq) + Cl-(aq)')
ChemicalEquation('NaCl(s) <=> Cl-(aq) + Na+(aq)')
>>> ChemicalSystem('E + S <=> ES -> E + P')
ChemicalSystem(["E + S <=> ES", "ES -> E + P"])

`atoms` is also available:

>>> from pyrrole import atoms
>>> acetate = atoms.read_cclib("data/acetate/acetate.out", "AcO-(g)")
"""

from pyrrole import atoms  # noqa
from pyrrole.core import ChemicalEquation, ChemicalSystem  # noqa
from pyrrole.drawing import draw_diagram  # noqa
