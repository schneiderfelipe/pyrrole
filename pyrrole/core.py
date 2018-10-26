#!/usr/bin/python3
# -*- encoding: utf-8 -*-

"""Tools for handling chemical equations and systems."""

import re
import numbers as _numbers
try:
    from collections.abc import Iterable as _Iterable
    from collections.abc import Mapping as _Mapping
    from collections.abc import Sequence as _Sequence
except ImportError:
    from collections import Iterable as _Iterable
    from collections import Mapping as _Mapping
    from collections import Sequence as _Sequence

import numpy as _np
import pandas as _pd
import networkx as _nx
import pyparsing as _pp


def _parse_chemical_equation(value):
    """
    Parse the chemical equation mini-language.

    See the docstring of `ChemicalEquation` for more.

    Parameters
    ----------
    value : `str`
        A string in chemical equation mini-language.

    Returns
    -------
    mapping
        A mapping in the format specified by the mini-language (see notes on
        `ChemicalEquation`).

    Examples
    --------
    >>> from pyrrole.core import _parse_chemical_equation
    >>> parsed = _parse_chemical_equation('4 A + 3 B <- 2 C + D')
    >>> parsed['arrow']
    '->'
    >>> parsed['products'][1]['species']
    'B'
    >>> parsed['reactants'][0]['coefficient']
    2

    """
    arrow = _pp.oneOf('-> <- <=>').setResultsName('arrow')
    species = _pp.Word(_pp.printables).setResultsName('species')
    coefficient = (_pp.Optional(_pp.Word(_pp.nums), default=1)
                   .setParseAction(_pp.tokenMap(int))
                   .setResultsName('coefficient'))
    group_ = _pp.Group(coefficient + _pp.Optional(_pp.Suppress('*')) + species)
    reactants = ((group_ + _pp.ZeroOrMore(_pp.Suppress('+') + group_))
                 .setResultsName('reactants'))
    products = ((group_ + _pp.ZeroOrMore(_pp.Suppress('+') + group_))
                .setResultsName('products'))

    grammar = reactants + arrow + products
    parsed = grammar.parseString(value).asDict()

    if parsed['arrow'] == '<-':
        parsed['reactants'], parsed['products'] \
            = parsed['products'], parsed['reactants']
        parsed['arrow'] = '->'

    return parsed


def _get_chemical_equation_piece(species_list, coefficients):
    """
    Produce a string from chemical species and their coefficients.

    Parameters
    ----------
    species_list : iterable of `str`
        Iterable of chemical species.
    coefficients : iterable of `float`
        Nonzero stoichiometric coefficients. The length of `species_list` and
        `coefficients` must be the same. Negative values are made positive and
        zeros are ignored along with their respective species.

    Examples
    --------
    >>> from pyrrole.core import _get_chemical_equation_piece
    >>> _get_chemical_equation_piece(["AcOH"], [2])
    '2 AcOH'
    >>> _get_chemical_equation_piece(["AcO-", "H+"], [-1, -1])
    'AcO- + H+'
    >>> _get_chemical_equation_piece("ABCD", [-2, -1, 0, -1])
    '2 A + B + D'

    """
    def _get_token(species, coefficient):
        if coefficient == 1:
            return '{}'.format(species)
        else:
            return '{:g} {}'.format(coefficient, species)

    bag = []
    for species, coefficient in zip(species_list, coefficients):
        if coefficient < 0:
            coefficient = -coefficient
        if coefficient > 0:
            bag.append(_get_token(species, coefficient))
    return '{}'.format(' + '.join(bag))


class ChemicalEquation:
    """
    An object for manipulating chemical equations in a way similar to vectors.

    This class provides an abstraction for chemical equations, generalizing
    equilibria, reactions, phase transitions and others. Conceptually,
    `ChemicalEquation` works like a vector that can be manipulated and operated
    upon. This allows the calculation of reaction free energies from data about
    chemical species, for instance.

    Parameters
    ----------
    value : `str`, `ChemicalEquation`, mapping or `Series`
        A string in chemical equation mini-language (see notes below), another
        `ChemicalEquation` object (which is copied) or either a mapping or
        `Series` from chemical species (`str`) to *signed* stoichiometric
        coefficients. See examples below.
    data : `pandas.DataFrame`, optional
        A `data` object, i.e., a table whose rows store information about
        chemical species, indexed by chemical species.
    arrow : `str`, optional
        Arrow symbol to use if `value` is a mapping or `Series`,
        ignored otherwise.

    Attributes
    ----------
    arrow : `str`
        Arrow symbol that separates `reactants` and `products` in
        a chemical equation. It is always equal to either ``'->'`` or
        ``'<=>'``.
    coefficient : mapping of `str` to `float`
        A mapping relating each `species` to its signed
        stoichiometric coefficient. Values are positive (negative) for
        `products` (`reactants`).
    products : iterable of `str`
        Chemical product species, i.e., those on the right-hand side of the
        equation.
    reactants : iterable of `str`
        Chemical reactant species, i.e., those on the left-hand side of the
        equation.
    species : iterable of `str`
        All species, i.e., union of all `reactants` and `products`.

    Notes
    -----
    Chemical equations in pyrrole are defined according to the following
    mini-language (white spaces are ignored)::

                   equation ::= reactants arrow products
        reactants, products ::=      coefficient ['*'] species
                                ['+' coefficient ['*'] species]*
                coefficient ::= [integers] (defaults to 1)
                    species ::= mix of printable characters
                      arrow ::= '->' | '<-' | '<=>'

    Examples
    --------
    >>> from pyrrole import ChemicalEquation
    >>> ChemicalEquation('H2O <=> H+ + OH-')
    ChemicalEquation('H2O <=> H+ + OH-')

    Chemical equations are stored in the usual order, even when given inverted
    (see notes above). This means that, although understood, ``'<-'`` is never
    a value of `arrow`:

    >>> equation = ChemicalEquation('2 CO <- CO2 + C')
    >>> equation.arrow
    '->'
    >>> equation  # stored in the usual order
    ChemicalEquation('C + CO2 -> 2 CO')

    Chemical species that appear in both reactants and products are simplified:

    >>> ChemicalEquation('A + B -> 2 A')
    ChemicalEquation('B -> A')

    Chemical equations can be manipulated as entities in a vector space. For
    instance, they can be added and subtracted, ...

    >>> (-ChemicalEquation('AgCl(s) <=> Ag+(aq) + Cl-(aq)')
    ...  +ChemicalEquation('NaCl(s) <=> Na+(aq) + Cl-(aq)')
    ...  +ChemicalEquation('AgNO3(s) <=> Ag+(aq) + NO3-(aq)')
    ...  -ChemicalEquation('NaNO3(s) <=> Na+(aq) + NO3-(aq)'))
    ChemicalEquation('AgNO3(s) + NaCl(s) <=> AgCl(s) + NaNO3(s)')

    ... and also divided and multiplied by numbers:

    >>> 10 * ChemicalEquation('CH4 + 2 O2 -> CO2 + 2 H2O') / 2
    ChemicalEquation('5 CH4 + 10 O2 -> 5 CO2 + 10 H2O')

    Stoichiometric coefficients are always *signed*, that is, positive for
    `products` and negative for `reactants`:

    >>> equation = ChemicalEquation('NaCl + AgNO3 -> NaNO3 + AgCl(v)')
    >>> equation.coefficient['AgCl(v)']  # a product
    1.0
    >>> equation.coefficient['AgNO3']  # a reactant
    -1.0

    Convenient attributes make it easy to obtain iterables of chemical species.

    >>> "AgCl(v)" in equation.reactants
    False
    >>> for species in equation.species:
    ...     print(species)
    AgCl(v)
    AgNO3
    NaCl
    NaNO3

    Although not recomended for everyday usage, a chemical equation can also be
    created from a complete set of stoichiometric coefficients. This makes it
    easy to convert other objects into a `ChemicalEquation` object:

    >>> ChemicalEquation({'NaCl': -1, 'AgNO3': -1,
    ...                   'NaNO3': 1, 'AgCl(v)': 1}, arrow='->')
    ChemicalEquation('AgNO3 + NaCl -> AgCl(v) + NaNO3')

    """

    def __init__(self, value, data=None, arrow=None):
        """See the docstring for this class."""
        # TODO: make tests for usage of data.
        if isinstance(value, ChemicalEquation):
            self.arrow = value.arrow
            self.coefficient = value.coefficient
            # TODO: make a test for this if.
            if data is None:
                data = value.data
        elif isinstance(value, str):
            parsed = _parse_chemical_equation(value)
            coefficient_products = _pd.Series({
                product['species']: +product['coefficient']
                for product in parsed['products']})
            coefficient_reactants = _pd.Series({
                reactant['species']: -reactant['coefficient']
                for reactant in parsed['reactants']})

            self.arrow = parsed['arrow']
            self.coefficient = coefficient_reactants.add(coefficient_products,
                                                         fill_value=0)
        elif isinstance(value, (_Mapping, _pd.Series)):
            # TODO: make test.
            if arrow not in {'->', '<=>'}:
                raise ValueError("arrow must be either '->' or '<=>' ('{}' "
                                 "given)".format(arrow))

            self.arrow = arrow
            self.coefficient = _pd.Series(value)
        else:
            raise TypeError("value must be either str, mapping or "
                            "Series ('{}' "
                            "given)".format(type(value).__name__))

        self.coefficient = self.coefficient.rename(self.__str__())
        self.data = data
        if self.data is not None and not isinstance(self.data, _pd.DataFrame):
            self.data = _pd.DataFrame(self.data)

    def _get_products(self):
        return self.coefficient[self.coefficient > 0].index
    products = property(_get_products)

    def _get_reactants(self):
        return self.coefficient[self.coefficient < 0].index
    reactants = property(_get_reactants)

    def _get_species(self):
        return self.coefficient.index
    species = property(_get_species)

    def __add__(self, other):
        """Add chemical equations as if they were vectors."""
        if not isinstance(other, ChemicalEquation):
            raise NotImplementedError

        return ChemicalEquation(self.coefficient.add(other.coefficient,
                                                     fill_value=0),
                                arrow=self.arrow)

    def __sub__(self, other):
        """Subtract chemical equations as if they were vectors."""
        if not isinstance(other, ChemicalEquation):
            raise NotImplementedError

        return ChemicalEquation(self.coefficient.sub(other.coefficient,
                                                     fill_value=0),
                                arrow=self.arrow)

    def __mul__(self, other):
        """Multiply chemical equation by a number."""
        if not isinstance(other, _numbers.Number):
            raise NotImplementedError

        return ChemicalEquation(self.coefficient.mul(other, fill_value=0),
                                arrow=self.arrow)
    __rmul__ = __mul__

    def __truediv__(self, other):
        """Divide chemical equation by a number."""
        if not isinstance(other, _numbers.Number):
            raise NotImplementedError

        return ChemicalEquation(self.coefficient.div(other, fill_value=0),
                                arrow=self.arrow)

    def __div__(self, other):
        """Ensure "true" division always takes place."""
        return self.__truediv__(other)

    def __pos__(self):
        """Make unary plus operator be equivalent to multiply by one."""
        return self

    def __neg__(self):
        """Make unary minus operator be equivalent to multiply by minus one."""
        return ChemicalEquation(self.coefficient.mul(-1, fill_value=0),
                                arrow=self.arrow)

    def __eq__(self, other):
        """Compare chemical equations in terms of coefficients."""
        diff = self.__sub__(other)
        return all(diff.coefficient == 0)

    def __repr__(self):
        """Build a string representation of this object."""
        return "ChemicalEquation('{}')".format(self)

    def __str__(self):
        """Build a unique, parsable string for this chemical equation."""
        products = []
        reactants = []
        for species, coefficient in sorted(self.coefficient.items()):
            if coefficient < 0:
                bag = reactants
            elif coefficient > 0:
                bag = products
            else:
                continue
            bag.append(_get_chemical_equation_piece([species], [coefficient]))

        return '{} {} {}'.format(' + '.join(reactants), self.arrow,
                                 ' + '.join(products))

    def to_series(self, only=None):
        """
        Produce a data record for `ChemicalEquation`.

        All possible linear differences for all numeric attributes are computed
        and stored in the returned `pandas.Series` object (see examples below).
        This allows for easy application and manipulation of
        `Hess's law <https://en.wikipedia.org/wiki/Hess%27s_law>`_ to chemical
        equations (see examples below).

        Parameters
        ----------
        only : ``"reactants"``, ``"products"``, optional
            Instead of the standard behaviour (difference of sums), sum numeric
            attributes of either reactants or products only. If given, absolute
            coefficients are used.

        Returns
        -------
        series : `pandas.Series`
            Data record of attribute differences, whose name is the canonical
            string representation of the `ChemicalEquation` or, if `only` is
            given, a string representing either reactants or products (see
            examples below).

        Examples
        --------
        >>> from pyrrole import ChemicalEquation
        >>> from pyrrole.atoms import create_data, read_cclib
        >>> data = create_data(
        ...     read_cclib("data/acetate/acetic_acid.out", "AcOH(g)"),
        ...     read_cclib("data/acetate/acetic_acid@water.out", "AcOH(aq)"))
        >>> equilibrium = ChemicalEquation("AcOH(g) <=> AcOH(aq)",
        ...                                data)
        >>> equilibrium.to_series()
        charge         0.000000
        enthalpy      -0.010958
        entropy       -0.000198
        freeenergy    -0.010759
        mult           0.000000
        natom          0.000000
        nbasis         0.000000
        nmo            0.000000
        pressure       0.000000
        temperature    0.000000
        Name: AcOH(g) <=> AcOH(aq), dtype: float64

        Sums of either reactants or products can be computed:

        >>> equilibrium.to_series("reactants")
        charge           0.000000
        enthalpy      -228.533374
        entropy          0.031135
        freeenergy    -228.564509
        mult             1.000000
        natom            8.000000
        nbasis          68.000000
        nmo             68.000000
        pressure         1.000000
        temperature    298.150000
        Name: AcOH(g), dtype: float64

        """
        if self.data is None:
            # TODO: should an empty Series be returned?
            raise ValueError("data not defined")

        # TODO: find a way to keep categorical columns. Keep if they match?
        columns = self.data.select_dtypes('number').columns

        if only is None:
            species = self.species
        elif only == "reactants":
            species = sorted(self.reactants)
        elif only == "products":
            species = sorted(self.products)
        else:
            raise ValueError("only must be either 'reactants' or 'products' "
                             "('{}' given)".format(only))

        if all([s in self.data.index for s in species]):
            # TODO: if two rows in self.data have the same indices, this must
            # sort them somehow (custom compare function?) and get the first
            # row
            series = (self.data.loc[species, columns]
                      .mul(self.coefficient, axis="index").sum("index"))
        else:
            series = _pd.Series(_np.nan, index=columns)

        if only is None:
            name = self.__str__()
        else:
            coefficients = self.coefficient[species]
            name = _get_chemical_equation_piece(species, coefficients)
            if only == "reactants":
                series = -series

        # Avoid negative zero
        # (see https://stackoverflow.com/a/11010791/4039050)
        series = series + 0.

        return series.rename(name)


def _split_arrows(value):
    """
    Split a string with sequential chemical equations into separate strings.

    Strings in odd positions in the returned iterable represent sums of
    chemical species (with possible stoichiometric coefficients). Strings in
    even positions represent arrow symbols. See examples below.

    Parameters
    ----------
    value : `str`
        A string with sequential chemical equations in the mini-language (see
        notes on `ChemicalEquation`).

    Returns
    -------
    iterable of `str`
        An iterable of strings. Odd positions represent sums of chemical
        species (with possible stoichiometric coefficients). Strings in even
        positions represent arrow symbols. See examples below.

    Notes
    -----
    Spaces are not striped from the returned strings (see examples below).

    Examples
    --------
    >>> from pyrrole.core import _split_arrows
    >>> _split_arrows('A -> B')
    ['A ', '->', ' B']

    """
    return re.split(r"(->|<-|<=>)", value)


def _split_chemical_equations(value):
    """
    Split a string with sequential chemical equations into separate strings.

    Each string in the returned iterable represents a single chemical equation
    of the input.
    See the docstrings of `ChemicalEquation` and `ChemicalSystem` for more.

    Parameters
    ----------
    value : `str`
        A string with sequential chemical equations in the mini-language (see
        notes on `ChemicalEquation`).

    Returns
    -------
    iterable of `str`
        An iterable of strings in the format specified by the mini-language
        (see notes on `ChemicalEquation`).

    Examples
    --------
    >>> from pyrrole.core import _split_chemical_equations
    >>> _split_chemical_equations('A + B -> C + D -> D + E <=> F + G <- H + I')
    ['A + B -> C + D', 'C + D -> D + E', 'D + E <=> F + G', 'F + G <- H + I']

    """
    pieces = _split_arrows(value)
    return [(pieces[i] +
             pieces[i + 1] +
             pieces[i + 2]).strip()
            for i in range(0, len(pieces) - 2, 2)]


# TODO: improve the examples for this class.
class ChemicalSystem:
    """
    Abstraction for models consisting of a set of chemical equations.

    Parameters
    ----------
    values : `ChemicalEquation`, `str`, sequence of `ChemicalEquation` or `str`
        Definitions of chemical equations. This can either be a single value or
        an iterable and accepts anything that can become a `ChemicalEquation`,
        or strings with consecutive equations (see examples below).
    data : `pandas.DataFrame`, optional
        A `data` object, i.e., a table whose rows store information about
        chemical species, indexed by chemical species.

    Attributes
    ----------
    equations : iterable of `ChemicalEquation`
        Stored `ChemicalEquation` objects.

    Examples
    --------
    >>> from pyrrole import ChemicalSystem
    >>> ChemicalSystem("A <=> B -> 2 C")
    ChemicalSystem(["A <=> B", "B -> 2 C"])

    Single chemical equations are also accepted, in which case the resulting
    model has a single equation:

    >>> ChemicalSystem(ChemicalEquation("A -> B"))
    ChemicalSystem(["A -> B"])

    Iterables can mix chemical equation definitions of different types:

    >>> ChemicalSystem(["A -> B", "A -> C <- D",
    ...                 ChemicalEquation("E -> A")])
    ChemicalSystem(["A -> B", "A -> C", "D -> C", "E -> A"])

    """

    def __init__(self, values, data=None):
        """See the docstring for this class."""
        if isinstance(values, str):
            self.equations = list(map(ChemicalEquation,
                                  _split_chemical_equations(values)))
        elif isinstance(values, (_Iterable, _Sequence)):
            self.equations = []
            for value in values:
                if isinstance(value, str):
                    self.equations.extend(map(ChemicalEquation,
                                          _split_chemical_equations(value)))
                else:
                    self.equations.append(ChemicalEquation(value))
        else:
            self.equations = [ChemicalEquation(values)]

        self.data = data
        if self.data is not None and not isinstance(self.data, _pd.DataFrame):
            self.data = _pd.DataFrame(self.data)
        for equation in self.equations:
            # TODO: make a test for this if.
            if equation.data is None:
                equation.data = self.data

    def __repr__(self):
        """Build a string representation of this object."""
        return "ChemicalSystem([{}])".format(
            ", ".join(['"' + str(equation) + '"'
                       for equation in self.equations]))

    def to_dataframe(self):
        """
        Produce a data table with records for all chemical equations.

        All possible differences for numeric attributes are computed and stored
        as columns in the returned `pandas.DataFrame` object (see examples
        below), whose rows represent chemical equations.

        In terms of behavior, this method can be seen as the `ChemicalEquation`
        counterpart of `create_data`.

        Returns
        -------
        dataframe : `pandas.DataFrame`
            Data table with records of attribute differences for every single
            `ChemicalEquation` object in the model.

        Examples
        --------
        >>> from pyrrole import ChemicalSystem
        >>> from pyrrole.atoms import create_data, read_cclib
        >>> data = create_data(
        ...     read_cclib("data/acetate/acetic_acid.out", "AcOH(g)"),
        ...     read_cclib("data/acetate/acetic_acid@water.out", "AcOH(aq)"))
        >>> data = data[["enthalpy", "freeenergy"]]
        >>> equilibrium = ChemicalSystem("AcOH(g) <=> AcOH(aq)", data)
        >>> equilibrium.to_dataframe()  # doctest: +NORMALIZE_WHITESPACE
                              enthalpy  freeenergy
        chemical_equation
        AcOH(g) <=> AcOH(aq) -0.010958   -0.010759

        """
        dataframe = _pd.DataFrame([equation.to_series()
                                   for equation in self.equations])
        dataframe.index.name = "chemical_equation"
        return dataframe

    def to_digraph(self):
        """
        Compute a directed graph for the chemical system.

        Returns
        -------
        digraph : `networkx.DiGraph`
            Graph nodes are reactants and/or products of chemical equations,
            while edges represent the equations themselves. Double ended edges
            are used to represent equilibria. Attributes are computed with
            `ChemicalEquation.to_series` for each equation (see examples
            below).

        Examples
        --------
        >>> from pyrrole import ChemicalSystem
        >>> from pyrrole.atoms import create_data, read_cclib
        >>> data = create_data(
        ...     read_cclib("data/acetate/acetic_acid.out", "AcOH(g)"),
        ...     read_cclib("data/acetate/acetic_acid@water.out", "AcOH(aq)"))
        >>> equilibrium = ChemicalSystem("AcOH(g) <=> AcOH(aq)", data)
        >>> digraph = equilibrium.to_digraph()
        >>> sorted(digraph.nodes(data='freeenergy'))
        [('AcOH(aq)', -228.57526805), ('AcOH(g)', -228.56450866)]
        >>> digraph.number_of_nodes()
        2
        >>> digraph.number_of_edges()
        2

        """
        # TODO: make test for this
        digraph = _nx.DiGraph()
        for equation in self.equations:
            reactants, arrow, products = [value.strip() for value
                                          in _split_arrows(str(equation))]

            try:
                attr = equation.to_series("reactants").to_dict()
            except ValueError:
                attr = dict()
            digraph.add_node(reactants, **attr)

            try:
                attr = equation.to_series("products").to_dict()
            except ValueError:
                attr = dict()
            digraph.add_node(products, **attr)

            try:
                attr = equation.to_series().to_dict()
            except ValueError:
                attr = dict()
            digraph.add_edge(reactants, products, **attr)
            if arrow == '<=>':
                digraph.add_edge(products, reactants, **attr)

        return digraph
