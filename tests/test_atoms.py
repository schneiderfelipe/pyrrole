#!/usr/bin/python3
# -*- encoding: utf-8 -*-

"""Tests for `pyrrole.atoms`."""

import numpy as np
import pandas as pd
from cclib import ccopen
from nose.tools import assert_equal
from numpy.testing import assert_allclose, assert_array_equal

from pyrrole.atoms import Atoms, create_data, read_cclib, read_pybel


def test_can_make_atoms_from_dict():
    """Test if we can make Atoms from dictionary."""
    atomcoords = [[0.0, 0.0, 0.0], [0.0, 0.0, 1.128]]
    atomnos = [6, 8]
    charge = 0

    mol = Atoms({
        'atomcoords': atomcoords,
        'atomnos': atomnos,
        'charge': charge
    })

    assert_allclose(mol.atomcoords[-1], atomcoords)
    assert_array_equal(mol.atomnos, atomnos)
    assert_equal(mol.charge, charge)


def test_can_make_atoms_from_series():
    """Test if we can make Atoms from Series."""
    atomcoords = [[0.0, 0.0, 0.0], [0.0, 0.0, 1.128]]
    atomnos = [6, 8]
    charge = 0

    attributes = {
        'atomcoords': atomcoords,
        'atomnos': atomnos,
        'charge': charge
    }

    mol1 = Atoms(attributes)
    mol2 = Atoms(pd.Series(attributes))

    assert_equal(mol1, mol2)


def test_can_make_atoms_with_read_cclib_and_ccopen():
    """Test if we can make Atoms with read_cclib and ccopen."""
    path = 'data/pyrrole.out'
    logfile = ccopen(path)

    mol = read_cclib(logfile.parse())

    assert_allclose(mol.atomcoords[-1],
                    [[-1.168043, 0.008417, 0.000141],
                     [-0.377180, 1.126730, 0.000400],
                     [0.970341, 0.688618, 0.000160],
                     [0.952121, -0.680756, -0.000226],
                     [-0.349506, -1.079198, -0.000239],
                     [-2.240180, -0.093370, 0.000197],
                     [-0.725206, 2.146546, 0.000720],
                     [1.851648, 1.308599, 0.000269],
                     [1.759012, -1.393872, -0.000499],
                     [-0.658645, -2.033234, -0.000482]])
    assert_allclose(mol.enthalpy, -209.61386780)
    assert_allclose(mol.entropy, 0.03005736)
    assert_allclose(mol.freeenergy, -209.64392516)
    assert_array_equal(mol.atomnos,
                       [6, 6, 6, 6, 7, 1, 1, 1, 1, 1])
    assert_equal(mol.charge, 0)
    assert_equal(mol.mult, 1)
    assert_equal(mol.natom, 10)

    # TODO: is this is the pure scf energy (no corrections added)?
    assert_equal(mol.scfenergies[-1], -5706.623393767155)


def test_can_access_attributes():
    """Test if we can access data as attributes."""
    path = 'data/pyrrole.out'
    mol = read_cclib(path)

    for key, value in mol.attributes.items():
        try:
            assert_equal(getattr(mol, key), value)
        except ValueError:
            assert_array_equal(getattr(mol, key), value)


def test_can_get_logfile_path():
    """Test if logfile path is obtainable."""
    path = 'data/pyrrole.out'
    mol = read_cclib(path)

    assert_equal(mol.jobfilename, path)


def test_does_read_cclib_call_parse():
    """Test if read_cclib calls parse on object if necessary."""
    path = 'data/pyrrole.out'
    logfile = ccopen(path)

    mol1 = read_cclib(logfile)
    mol2 = read_cclib(logfile.parse())

    # Make a fair comparison.
    mol1.attributes['name'] = None
    mol1.attributes['jobfilename'] = None

    assert_equal(mol1, mol2)


def test_can_make_atoms_with_read_cclib_and_path():
    """Test if we can make Atoms with read_cclib and path to logfile."""
    path = 'data/pyrrole.out'
    logfile = ccopen(path)

    mol1 = read_cclib(path)
    mol2 = read_cclib(logfile)

    assert_equal(mol1, mol2)


def test_can_print_atoms():
    """Test if we can print Atom instances."""
    atomcoords = [[0.0, 0.0, 0.0], [0.0, 0.0, 1.128]]
    atomnos = [6, 8]
    charge = 0
    mult = 1

    mol = Atoms({
        'atomcoords': atomcoords,
        'atomnos': atomnos,
        'charge': charge,
        'mult': mult
    })

    assert_equal(str(mol), """C          0.00000        0.00000        0.00000
O          0.00000        0.00000        1.12800""")


def test_can_print_atoms_with_scfenergy():
    """Test if we can print Atom instances with scfenergy."""
    path = 'data/pyrrole.out'
    mol = read_cclib(path)
    scfenergy = mol.scfenergies[-1]
    natom = mol.natom

    assert_equal(mol.to_string("xyz"), str(mol))
    assert_equal(mol.to_string("xyz", with_header=True), """{}
{}, scfenergy={} eV
C         -1.16804        0.00842        0.00014
C         -0.37718        1.12673        0.00040
C          0.97034        0.68862        0.00016
C          0.95212       -0.68076       -0.00023
N         -0.34951       -1.07920       -0.00024
H         -2.24018       -0.09337        0.00020
H         -0.72521        2.14655        0.00072
H          1.85165        1.30860        0.00027
H          1.75901       -1.39387       -0.00050
H         -0.65865       -2.03323       -0.00048""".format(natom, path,
                                                           scfenergy))


def test_can_print_atoms_in_formats_and_dialects():
    """Test if we can write Atom instances in many formats and dialects."""
    water_dimer = read_pybel("data/water-dimer.xyz")

    assert_equal(water_dimer.to_string("xyz"),
                 water_dimer.to_string("xyz", dialect="standard"))
    assert_equal(water_dimer.to_string("xyz", with_header=True),
                 """6
data/water-dimer.xyz
O         -1.62893       -0.04138        0.37137
H         -0.69803       -0.09168        0.09337
H         -2.06663       -0.73498       -0.13663
O          1.21457        0.03172       -0.27623
H          1.44927        0.91672       -0.58573
H          1.72977       -0.08038        0.53387""")

    assert_equal(water_dimer.to_string("xyz"),
                 """O         -1.62893       -0.04138        0.37137
H         -0.69803       -0.09168        0.09337
H         -2.06663       -0.73498       -0.13663
O          1.21457        0.03172       -0.27623
H          1.44927        0.91672       -0.58573
H          1.72977       -0.08038        0.53387""")

    assert_equal(water_dimer.to_string("xyz"),
                 water_dimer.to_string("xyz", dialect="adf"))
    assert_equal(water_dimer.to_string("xyz", dialect="ADF",
                                       fragment_id="water"),
                 """O         -1.62893       -0.04138        0.37137       f=water
H         -0.69803       -0.09168        0.09337       f=water
H         -2.06663       -0.73498       -0.13663       f=water
O          1.21457        0.03172       -0.27623       f=water
H          1.44927        0.91672       -0.58573       f=water
H          1.72977       -0.08038        0.53387       f=water""")

    assert_equal(water_dimer.to_string("gamin"),
                 """O      8.0     -1.6289300000   -0.0413800000    0.3713700000
H      1.0     -0.6980300000   -0.0916800000    0.0933700000
H      1.0     -2.0666300000   -0.7349800000   -0.1366300000
O      8.0      1.2145700000    0.0317200000   -0.2762300000
H      1.0      1.4492700000    0.9167200000   -0.5857300000
H      1.0      1.7297700000   -0.0803800000    0.5338700000""")

    assert_equal(water_dimer.to_string("mop", with_header=True),
                 """PUT KEYWORDS HERE
data/water-dimer.xyz
O  -1.62893 1 -0.04138 1  0.37137 1
H  -0.69803 1 -0.09168 1  0.09337 1
H  -2.06663 1 -0.73498 1 -0.13663 1
O   1.21457 1  0.03172 1 -0.27623 1
H   1.44927 1  0.91672 1 -0.58573 1
H   1.72977 1 -0.08038 1  0.53387 1""")

    assert_equal(water_dimer.to_string("mop", constraints=[i for i in range(4)
                                                           if i % 2 == 0]),
                 """O  -1.62893 0 -0.04138 0  0.37137 0
H  -0.69803 1 -0.09168 1  0.09337 1
H  -2.06663 0 -0.73498 0 -0.13663 0
O   1.21457 1  0.03172 1 -0.27623 1
H   1.44927 1  0.91672 1 -0.58573 1
H   1.72977 1 -0.08038 1  0.53387 1""")


def test_can_convert_atoms_to_series():
    """Test if we can convert Atom instances to series."""
    atomcoords = [[0., 0., 0.], [0., 0., 1.21]]
    atomnos = [8, 8]
    charge = 0
    mult = 3
    name = 'dioxygen'

    mol = Atoms({
        'atomcoords': atomcoords,
        'atomnos': atomnos,
        'charge': charge,
        'mult': mult,
        'name': name
    })

    series = mol.to_series()
    assert_allclose(series['atomcoords'][-1], atomcoords)
    assert_array_equal(series['atomnos'], atomnos)
    assert_equal(series['charge'], charge)
    assert_equal(series['mult'], mult)
    assert_equal(series['name'], name)


def test_can_create_data_from_atoms():
    """Test if we can create a data object from Atom instances."""
    pyrrole_name = 'pyrrole'
    pyrrole_path = 'data/pyrrole.out'
    pyrrole_mol = read_cclib(pyrrole_path, pyrrole_name)

    pyrrolate_path = 'data/pyrrolate.out'
    pyrrolate_mol = read_cclib(pyrrolate_path)

    empty_mol = Atoms({})

    dioxygen_name = 'dioxygen'
    dioxygen_atomcoords = [[0., 0., 0.], [0., 0., 1.21]]
    dioxygen_atomnos = [8, 8]
    dioxygen_mult = 3
    dioxygen_mol = Atoms({
        'name': dioxygen_name,
        'atomcoords': dioxygen_atomcoords,
        'atomnos': dioxygen_atomnos,
        'mult': dioxygen_mult,
    })

    cyanide_name = 'cyanide'
    cyanide_atomcoords = [[0., 0., 0.], [0., 0., 1.136]]
    cyanide_atomnos = [7, 8]
    cyanide_charge = -1
    cyanide_mol = Atoms({
        'name': cyanide_name,
        'atomcoords': cyanide_atomcoords,
        'atomnos': cyanide_atomnos,
        'charge': cyanide_charge,
    })

    data = create_data([pyrrole_mol,
                        pyrrolate_mol,
                        empty_mol,
                        dioxygen_mol,
                        cyanide_mol])

    assert_equal(data.loc[pyrrole_name, 'jobfilename'], pyrrole_path)
    assert_array_equal(data.loc[pyrrole_name, 'atomnos'],
                       [6, 6, 6, 6, 7, 1, 1, 1, 1, 1])
    assert_equal(data.loc[pyrrole_name, 'charge'], 0.)
    assert_equal(data.loc[pyrrole_name, 'mult'], 1.)

    assert_equal(data.loc[pyrrolate_path, 'jobfilename'], pyrrolate_path)
    assert_array_equal(data.loc[pyrrolate_path, 'atomnos'],
                       [6, 6, 6, 6, 7, 1, 1, 1, 1])
    assert_equal(data.loc[pyrrolate_path, 'charge'], -1.)
    assert_equal(data.loc[pyrrolate_path, 'mult'], 1.)

    assert(np.all(pd.isna(data.loc[np.nan])))

    assert(pd.isna(data.loc[dioxygen_name, 'jobfilename']))
    assert_allclose(data.loc[dioxygen_name, 'atomcoords'][-1],
                    dioxygen_atomcoords)
    assert_array_equal(data.loc[dioxygen_name, 'atomnos'], dioxygen_atomnos)
    assert(pd.isna(data.loc[dioxygen_name, 'charge']))
    assert_equal(data.loc[dioxygen_name, 'mult'], dioxygen_mult)

    assert(pd.isna(data.loc[cyanide_name, 'jobfilename']))
    assert_allclose(data.loc[cyanide_name, 'atomcoords'][-1],
                    cyanide_atomcoords)
    assert_array_equal(data.loc[cyanide_name, 'atomnos'], cyanide_atomnos)
    assert_equal(data.loc[cyanide_name, 'charge'], cyanide_charge)
    assert(pd.isna(data.loc[cyanide_name, 'mult']))


def test_can_fragment_molecule():
    """Test if we can fragment a molecule."""
    water_dimer = read_pybel("data/water-dimer.xyz")

    assert_equal([frag.to_string("xyz", with_header=True)
                  for frag in water_dimer.split()],
                 ["""3

O         -1.62893       -0.04138        0.37137
H         -0.69803       -0.09168        0.09337
H         -2.06663       -0.73498       -0.13663""",
                  """3

O          1.21457        0.03172       -0.27623
H          1.72977       -0.08038        0.53387
H          1.44927        0.91672       -0.58573"""])

    assert_equal([frag.to_string("xyz", with_header=True)
                  for frag in water_dimer.split([range(6)])],
                 ["""6

O         -1.62893       -0.04138        0.37137
H         -0.69803       -0.09168        0.09337
H         -2.06663       -0.73498       -0.13663
O          1.21457        0.03172       -0.27623
H          1.44927        0.91672       -0.58573
H          1.72977       -0.08038        0.53387"""])

    assert_equal([frag.to_string("xyz", with_header=True)
                  for frag in water_dimer.split([(2, 1), [5]])],
                 ["""2

H         -2.06663       -0.73498       -0.13663
H         -0.69803       -0.09168        0.09337""",
                  """1

H          1.72977       -0.08038        0.53387"""])


