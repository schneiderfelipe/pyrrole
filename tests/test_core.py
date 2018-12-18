#!/usr/bin/python3
# -*- encoding: utf-8 -*-

"""Tests for `pyrrole.core`."""

import numpy as np
import pandas as pd
from nose.tools import assert_equal
from numpy.testing import assert_allclose

from pyrrole.core import (_parse_chemical_equation, _split_arrows,
                          ChemicalEquation, ChemicalSystem)
from pyrrole.atoms import create_data, read_cclib


def test_can_parse_chemical_equation_the_low_level_way():
    """Test if we can parse chemical equations."""
    parsed = _parse_chemical_equation('4 A(aq) + 3 B <- 2 C* + D#')
    expected = {'reactants': [{'coefficient': 2, 'species': 'C*'},
                              {'coefficient': 1, 'species': 'D#'}],
                'products': [{'coefficient': 4, 'species': 'A(aq)'},
                             {'coefficient': 3, 'species': 'B'}],
                'arrow': '->'}
    assert_equal(parsed, expected)


def test_can_split_arrows_chemical_equation_the_low_level_way():
    """Test if we can split arrows in chemical equations."""
    parsed = _split_arrows('A + B -> C + D -> D + E <=> F + G <- H + I')
    expected = ['A + B ', '->', ' C + D ', '->', ' D + E ', '<=>', ' F + G ',
                '<-', ' H + I']
    assert_equal(parsed, expected)


def test_can_make_chemical_equation_from_str():
    """Test if we can create chemical equation from string."""
    equation = ChemicalEquation('4 A(aq) + 3 B <- 2 C* + D#')
    assert_equal(equation.arrow, '->')
    assert_equal(equation.coefficient['A(aq)'], 4)
    assert_equal(equation.coefficient['C*'], -2)
    assert_equal(set(equation.products), {'A(aq)', 'B'})
    assert_equal(set(equation.reactants), {'C*', 'D#'})
    assert_equal(set(equation.species), {'A(aq)', 'B', 'C*', 'D#'})


def test_can_make_chemical_equation_from_dict():
    """Test if we can create chemical equation from dictionary."""
    equation1 = ChemicalEquation({'A(aq)': 4, 'B': 3, 'C*': -2,
                                  'D#': -1}, arrow='->')
    equation2 = ChemicalEquation('4 A(aq) + 3 B <- 2 C* + D#')
    assert_equal(equation1, equation2)


def test_can_make_chemical_equation_from_series():
    """Test if we can create chemical equation from pandas.Series."""
    equation1 = ChemicalEquation(pd.Series({'A(aq)': 4, 'B': 3,
                                            'C*': -2, 'D#': -1}),
                                 arrow='->')
    equation2 = ChemicalEquation('4 A(aq) + 3 B <- 2 C* + D#')
    assert_equal(equation1, equation2)


def test_can_manipulate_chemical_equations_as_vectors():
    """Test if chemical equation can be manipulated as vectors."""
    equation1 = ChemicalEquation('A -> 2 B') / 2
    equation2 = ChemicalEquation('2 B -> C')
    equation3 = 2 * equation1 + equation2
    equation4 = -equation2 - (2 * equation1)
    assert_equal(equation3, ChemicalEquation('A -> C'))
    assert_equal(equation4, ChemicalEquation('C -> A'))


def test_can_parse_printed_chemical_equation():
    """Test if printed chemical equations can be correctly parsed."""
    equation1 = ChemicalEquation('4 A(aq) + 3 B <- 2 C* + D#')
    assert_equal(str(equation1), '2 C* + D# -> 4 A(aq) + 3 B')

    equation2 = ChemicalEquation(str(equation1))
    assert_equal(equation1, equation2)


def test_can_convert_chemical_equation_to_series():
    """Test if we can convert chemical equation to pandas.Series."""
    data = create_data(read_cclib("data/acetate/acetic_acid.out", "AcOH(g)"),
                       read_cclib("data/acetate/acetic_acid@water.out",
                                  "AcOH(aq)"))
    equilibrium = ChemicalEquation("AcOH(g) <=> AcOH(aq)", data)
    series = equilibrium.to_series()
    assert_equal(series['natom'], 0.)
    assert_equal(series['charge'], 0.)
    assert_equal(series['mult'], 0.)
    assert_equal(series['temperature'], 298.15)
    assert_equal(series['pressure'], 1.)
    assert_equal(series['nmo'], 0.)
    assert_equal(series['nbasis'], 0.)
    assert_allclose(series['enthalpy'], -0.01095788)
    assert_allclose(series['entropy'], -0.00019848)
    assert_allclose(series['freeenergy'], -0.01075939)


def test_is_chemical_equation_to_series_consistent():
    """Test if conversion of ChemicalEquation to pandas.Series is consistent."""  # noqa
    data = create_data(read_cclib("data/acetate/acetic_acid.out", "AcOH(g)"),
                       read_cclib("data/acetate/acetic_acid@water.out",
                                  "AcOH(aq)"))
    equilibrium = ChemicalEquation("AcOH(g) <=> AcOH(aq)", data)
    assert_allclose(equilibrium.to_series(intensive_columns=[]),
                    equilibrium.to_series("products")
                    - equilibrium.to_series("reactants"))


def test_can_obtain_species_from_chemical_equation():
    """Test if we can obtain species from chemical equations."""
    equation = ChemicalEquation('NaCl + AgNO3 -> NaNO3 + AgCl(v)')
    assert_equal(sorted(equation.species),
                 [u'AgCl(v)', u'AgNO3', u'NaCl', u'NaNO3'])
    assert_equal(sorted(equation.products),
                 [u'AgCl(v)', u'NaNO3'])
    assert_equal(sorted(equation.reactants),
                 [u'AgNO3', u'NaCl'])


def test_can_convert_chemical_system_to_dataframe():
    """Test if we can convert chemical system to dataframe."""
    data = create_data(read_cclib("data/acetate/acetic_acid.out", "AcOH(g)"),
                       read_cclib("data/acetate/acetic_acid@water.out",
                                  "AcOH(aq)"))
    equilibrium = ChemicalSystem("AcOH(g) <=> AcOH(aq)", data)
    dataframe = equilibrium.to_dataframe()
    assert(np.all(dataframe['natom'] == 0.))
    assert(np.all(dataframe['charge'] == 0.))
    assert(np.all(dataframe['mult'] == 0.))
    assert(np.all(dataframe['temperature'] == 298.15))
    assert(np.all(dataframe['pressure'] == 1.))
    assert(np.all(dataframe['nmo'] == 0.))
    assert(np.all(dataframe['nbasis'] == 0.))
    assert_allclose(dataframe['enthalpy'], -0.01095788)
    assert_allclose(dataframe['entropy'], -0.00019848)
    assert_allclose(dataframe['freeenergy'], -0.01075939)
