#!/usr/bin/python3
# -*- encoding: utf-8 -*-

"""Tests for `pyrrole.drawing`."""

import numpy as np
import pandas as pd
from nose.tools import assert_equal
from numpy.testing import assert_allclose

from pyrrole import ChemicalSystem
from pyrrole.atoms import create_data, read_cclib
from pyrrole.drawing import diagram_layout, tower_layout


def test_can_produce_tower_layout():
    """Test if we can produce a tower layout."""
    data = create_data(read_cclib("data/acetate/acetic_acid.out", "AcOH(g)"),
                       read_cclib("data/acetate/acetic_acid@water.out",
                                  "AcOH(aq)"))
    digraph = (ChemicalSystem("AcOH(g) <=> AcOH(aq)", data).to_digraph())

    layout = tower_layout(digraph)
    expected = {'AcOH(aq)': np.array([0., -228.57526805]),
                'AcOH(g)': np.array([0., -228.56450866])}

    assert_equal(expected.keys(), layout.keys())
    for name in expected:
        assert_allclose(layout[name], expected[name])


def test_can_produce_tower_layout_with_scale():
    """Test if we can produce a tower layout with scale."""
    data = create_data(read_cclib("data/acetate/acetic_acid.out", "AcOH(g)"),
                       read_cclib("data/acetate/acetic_acid@water.out",
                                  "AcOH(aq)"))
    digraph = (ChemicalSystem("AcOH(g) <=> AcOH(aq)", data).to_digraph())

    layout = tower_layout(digraph, scale=1)
    expected = {'AcOH(g)': np.array([0.,  1.]),
                'AcOH(aq)': np.array([0., -1.])}

    assert_equal(expected.keys(), layout.keys())
    for name in expected:
        assert_allclose(layout[name], expected[name])


def test_can_produce_diagram_layout():
    """Test if we can produce a diagram layout."""
    data = (pd.DataFrame([{"name": "Separated_Reactants", "freeenergy": 0.},
                          {"name": "mlC1", "freeenergy": -5.4},
                          {"name": "mlC2", "freeenergy": -15.6},
                          {"name": "mTS1", "freeenergy": 28.5, "color": "g"},
                          {"name": "mCARB1", "freeenergy": -9.7},
                          {"name": "mCARB2", "freeenergy": -19.8},
                          {"name": "mCARBX", "freeenergy": 20}])
            .set_index("name"))
    system = ChemicalSystem(
        ["Separated_Reactants -> mlC1 -> mTS1",
         "Separated_Reactants -> mlC2 -> mTS1",
         "mCARB2 <- mTS1 -> mCARB1",
         "Separated_Reactants -> mCARBX"], data)
    digraph = system.to_digraph()

    layout = diagram_layout(digraph)
    expected = {'Separated_Reactants': np.array([0., 0.]),
                'mlC1': np.array([1., -5.4]),
                'mTS1': np.array([2.,  28.5]),
                'mlC2': np.array([1., -15.6]),
                'mCARB2': np.array([3., -19.8]),
                'mCARB1': np.array([3., -9.7]),
                'mCARBX': np.array([1.,  20.])}

    assert_equal(expected.keys(), layout.keys())
    for name in expected:
        assert_allclose(layout[name], expected[name])


def test_can_produce_diagram_layout_with_scale():
    """Test if we can produce a diagram layout with scale."""
    data = (pd.DataFrame([{"name": "Separated_Reactants", "freeenergy": 0.},
                          {"name": "mlC1", "freeenergy": -5.4},
                          {"name": "mlC2", "freeenergy": -15.6},
                          {"name": "mTS1", "freeenergy": 28.5, "color": "g"},
                          {"name": "mCARB1", "freeenergy": -9.7},
                          {"name": "mCARB2", "freeenergy": -19.8},
                          {"name": "mCARBX", "freeenergy": 20}])
            .set_index("name"))
    system = ChemicalSystem(
        ["Separated_Reactants -> mlC1 -> mTS1",
         "Separated_Reactants -> mlC2 -> mTS1",
         "mCARB2 <- mTS1 -> mCARB1",
         "Separated_Reactants -> mCARBX"], data)
    digraph = system.to_digraph()

    layout = diagram_layout(digraph, scale=1)
    expected = {'Separated_Reactants': np.array([-0.05459057,  0.00992556]),
                'mlC1': np.array([-0.01985112, -0.17766749]),
                'mTS1': np.array([0.01488834,  1.]),
                'mlC2': np.array([-0.01985112, -0.53200993]),
                'mCARB2': np.array([0.04962779, -0.67791563]),
                'mCARB1': np.array([0.04962779, -0.32704715]),
                'mCARBX': np.array([-0.01985112,  0.70471464])}

    assert_equal(expected.keys(), layout.keys())
    for name in expected:
        assert_allclose(layout[name], expected[name], rtol=1e-06)
