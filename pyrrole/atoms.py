#!/usr/bin/python3
# -*- encoding: utf-8 -*-

"""Tools for handling molecular structures and their attributes."""

import os as _os
import re as _re

import numpy as _np
import pandas as _pd
import cclib as _cclib
import pybel as _pb
import openbabel as _ob
from cclib.parser import data as _data
from cclib.parser import logfileparser as _logfileparser
from cclib.bridge.cclib2openbabel import makecclib as _makecclib
from cclib.bridge.cclib2openbabel import makeopenbabel as _makeopenbabel


class Atoms:
    """
    Abstraction of fixed molecular structures as sequences of atoms.

    Properties are stored in `Atoms.attributes`, which follow naming and unit
    conventions of the `cclib project <https://cclib.github.io/>`_ as close as
    possible (see `read_cclib` for producing `Atoms` out of cclib objects or
    logfiles).

    Parameters
    ----------
    attributes : mapping or `pandas.Series`
        A mapping containing data attributes whose keys follow the naming
        conventions of the `cclib project <https://cclib.github.io/>`_ for
        properties as close as possible. Any mapping of attributes that can be
        converted to a `pandas.Series` object is accepted.

    Attributes
    ----------
    attribute : `pandas.Series`
        Molecular properties. Keys follow the naming practice of the `cclib
        project <https://cclib.github.io/>`_ for properties as close as
        possible. The object is named during initialization after `Atoms.name`
        (see notes below).

    Notes
    -----
    This class is intended to be used directly in very simple cases only. For
    complex situations, use factory functions such as `read_cclib`.

    `Atoms.name` is set to `Atoms.jobfilename` if `Atoms.name` does not exist.

    Examples
    --------
    An `Atoms` object can be instantiated by giving any mapping of
    attributes that can be converted to a `pandas.Series` object:

    >>> from pyrrole.atoms import Atoms
    >>> dioxygen = Atoms({
    ...     'atomcoords': [[0., 0., 0.],
    ...                    [0., 0., 1.21]],
    ...     'charge': 0,
    ...     'mult': 3,
    ...     'atomnos': [8, 8]})

    Attributes can be accessed by their keys, as given during initialization:

    >>> dioxygen.mult
    3
    >>> dioxygen.atomcoords[-1][1, 2]
    1.21

    Two `Atoms` objects are considered equal if their attributes match:

    >>> dioxygen_copy = Atoms({
    ...     'atomcoords': [[0., 0., 0.],
    ...                    [0., 0., 1.21]],
    ...     'charge': 0,
    ...     'mult': 3,
    ...     'atomnos': [8, 8]})
    >>> dioxygen == dioxygen_copy
    True
    >>> dioxygen is dioxygen_copy
    False

    If printed, `Atoms` objects produce coordinates in
    `xyz format <https://en.wikipedia.org/wiki/XYZ_file_format>`_:

    >>> from pyrrole.atoms import read_pybel
    >>> print(read_pybel("data/pyrrole.xyz", "pyrrole"))
    C         -1.11468        0.00000        0.33059
    C         -0.70848        0.00000       -0.97739
    C          0.70848        0.00000       -0.97739
    C          1.11468        0.00000        0.33059
    N          0.00000        0.00000        1.11189
    H         -2.10267        0.00000        0.75908
    H         -1.35484        0.00000       -1.83956
    H          1.35484        0.00000       -1.83956
    H          2.10267        0.00000        0.75908
    H          0.00000        0.00000        2.11476
    >>> print(dioxygen)
    O          0.00000        0.00000        0.00000
    O          0.00000        0.00000        1.21000

    """

    def __init__(self, attributes):
        """See the docstring for this class."""
        # TODO: make idempotent (allow receiving Atoms objects).
        self.attributes = _pd.Series(attributes)

        # TODO: properly document used attributes.
        if "name" in self.attributes:
            # TODO: make a test for this if.
            self.attributes.rename(self.name, inplace=True)
        elif "jobfilename" in self.attributes:
            self.attributes.rename(self.jobfilename, inplace=True)

        if "atomcoords" in self.attributes:
            self.attributes["atomcoords"] = _np.asanyarray(self.atomcoords)
            if len(self.atomcoords.shape) < 3:
                self.attributes["atomcoords"] = [self.atomcoords]

        self.name = self.attributes.name

    def __getattr__(self, value):
        """Wrap `Atoms.value` into `Atoms.attributes['value']`."""
        return self.attributes[value]

    def __eq__(self, other):
        """Compare molecules in terms of attributes."""
        # Expensive but definitely correct.
        return self.attributes.to_json() == other.attributes.to_json()

    def __str__(self):
        """Build a string with Cartesian coordinates for this molecule."""
        # TODO: make a read_string that is able to read the results of this
        # function.
        return self.to_string('xyz')

    def split(self, pattern=None):
        r"""
        Break molecule up into constituent fragments.

        By default (i.e., if `pattern` is `None`), each disconnected fragment
        is returned as a separate new `Atoms` object. This uses OpenBabel
        (through `OBMol.Separate`) and might not preserve atom order, depending
        on your version of the library.

        Parameters
        ----------
        pattern : iterable of iterable of `int`, optional
            Groupings of atoms into molecule fragments. Each element of
            `pattern` should be an iterable whose members are atom indices (see
            example below).

        Returns
        -------
        fragments : iterable of `Atoms`

        Examples
        --------
        >>> from pyrrole import atoms
        >>> water_dimer = atoms.read_pybel("data/water-dimer.xyz")

        "Natural fragmentation" is the default behaviour, i.e. all disconnected
        fragments are returned:

        >>> for frag in water_dimer.split():
        ...     print("{}\n".format(frag))
        O         -1.62893       -0.04138        0.37137
        H         -0.69803       -0.09168        0.09337
        H         -2.06663       -0.73498       -0.13663
        <BLANKLINE>
        O          1.21457        0.03172       -0.27623
        H          1.72977       -0.08038        0.53387
        H          1.44927        0.91672       -0.58573
        <BLANKLINE>

        Precise fragment grouping can be achieved by explicitly indicating
        which atoms belong to which fragments:

        >>> for frag in water_dimer.split([range(3), (5, 4), [3]]):
        ...     print("{}\n".format(frag))
        O         -1.62893       -0.04138        0.37137
        H         -0.69803       -0.09168        0.09337
        H         -2.06663       -0.73498       -0.13663
        <BLANKLINE>
        H          1.72977       -0.08038        0.53387
        H          1.44927        0.91672       -0.58573
        <BLANKLINE>
        O          1.21457        0.03172       -0.27623
        <BLANKLINE>

        """
        molecule_pybel = self.to_pybel()

        if pattern is None:
            fragments = [read_pybel(frag)
                         for frag in molecule_pybel.OBMol.Separate()]
        else:
            fragments = []
            for group in pattern:
                fragment_obmol = _pb.ob.OBMol()
                for i in group:
                    obatom = molecule_pybel.OBMol.GetAtomById(i)
                    fragment_obmol.InsertAtom(obatom)

                fragments.append(fragment_obmol)

            fragments = [read_pybel(frag) for frag in fragments]

        return fragments

    def to_series(self):
        """
        Produce a data record of attributes as a `pandas.Series` object.

        This method is useful for producing many data tables in pyrrole. See
        `ChemicalEquation.to_series` and `ChemicalSystem.to_dataframe` for
        more.

        Returns
        -------
        series : `pandas.Series`
            Data record of attributes.

        Examples
        --------
        >>> from pyrrole.atoms import Atoms
        >>> dioxygen = Atoms({'atomcoords': [[0., 0., 0.],
        ...                                  [0., 0., 1.21]],
        ...                   'atomnos': [8, 8],
        ...                   'mult': 3,
        ...                   'name': 'dioxygen'})
        >>> dioxygen.to_series()
        atomcoords    [[[0.0, 0.0, 0.0], [0.0, 0.0, 1.21]]]
        atomnos                                      [8, 8]
        mult                                              3
        name                                       dioxygen
        Name: dioxygen, dtype: object

        This method can be used to produce a copy of an `Atoms` object:

        >>> dioxygen_copy = Atoms(dioxygen.to_series())
        >>> dioxygen == dioxygen_copy
        True
        >>> dioxygen is dioxygen_copy
        False

        """
        return self.attributes

    def to_pybel(self):
        """
        Produce a Pybel Molecule object.

        It is based on the capabilities of OpenBabel through Pybel. The present
        object must have at least `atomcoords`, `atomnos`, `charge` and `mult`
        defined.

        Returns
        -------
        `pybel.Molecule`

        Examples
        --------
        >>> from pyrrole.atoms import Atoms
        >>> dioxygen = Atoms({'atomcoords': [[0., 0., 0.],
        ...                                  [0., 0., 1.21]],
        ...                   'atomnos': [8, 8],
        ...                   'charge': 0,
        ...                   'mult': 3,
        ...                   'name': 'dioxygen'})
        >>> mol = dioxygen.to_pybel()
        >>> mol.molwt
        31.9988

        """
        # TODO: This only exports last geometry by default.
        obmol = _makeopenbabel(self.atomcoords, self.atomnos, self.charge,
                               self.mult)

        title = self.name or ""
        if 'scfenergies' in self.attributes:
            title += ", scfenergy={} eV".format(self.scfenergies[-1])
        obmol.SetTitle(title)

        # TODO: make a test for this function.
        return _pb.Molecule(obmol)

    def to_string(self, format="smi", dialect=None, with_header=False,
                  fragment_id=None, constraints=None):
        r"""
        Produce a string representation of the molecule.

        This function wraps and extends the functionality of OpenBabel (which
        is accessible through `to_pybel`). Many chemical formats can thus be
        output (see the `pybel.outformats` variable for a list of available
        output formats).

        Parameters
        ----------
        format : `str`, optional
            Chemical file format of the returned string representation (see
            examples below).
        dialect : `str`, optional
            Format dialect. This encompasses enhancements provided for some
            subformats. If ``"standard"`` or `None`, the output provided by
            OpenBabel is used with no or minimal modification. See notes below.
        with_header : `bool`, optional
            If `format` encompasses a header, allow it in the returned string.
            This would be, for instance, the first two lines of data for
            ``format="xyz"`` (see examples below). This might not work with all
            dialects and/or formats.
        fragment_id : `str`, optional
            Indentify molecular fragments (see examples below). This might not
            work with all dialects and/or formats.
        constraints : iterable object of `int`
            Set cartesian constraints for selected atoms (see examples below).
            This might not work with all dialects and/or formats.

        Returns
        -------
        `str`
            String representation of molecule in the specified format and/or
            dialect.

        Notes
        -----
        Format dialects are subformats that support extended functionality.
        Currently supported dialects are:

        - for ``format="xyz"``:
            - ``"ADF"``, ``"ORCA"``.

        Examples
        --------
        >>> from pyrrole import atoms
        >>> dioxygen = atoms.Atoms({'atomcoords': [[0., 0., 0.],
        ...                                        [0., 0., 1.21]],
        ...                         'atomnos': [8, 8],
        ...                         'charge': 0,
        ...                         'mult': 3,
        ...                         'name': 'dioxygen'})

        By default, a SMILES string is returned:

        >>> dioxygen.to_string()
        'O=O\tdioxygen'

        Cartesian coordinates can be produced with ``format="xyz"``, which is
        equivalent to printing an `Atoms` instance:

        >>> print(dioxygen.to_string("xyz"))
        O          0.00000        0.00000        0.00000
        O          0.00000        0.00000        1.21000
        >>> print(dioxygen)
        O          0.00000        0.00000        0.00000
        O          0.00000        0.00000        1.21000

        Header lines are disabled by default (for ``format="xyz"``, for
        example, the header stores the number of atoms in the molecule and a
        comment or title line), but this can be reversed with
        ``with_header=True``:

        >>> print(dioxygen.to_string("xyz", with_header=True))
        2
        dioxygen
        O          0.00000        0.00000        0.00000
        O          0.00000        0.00000        1.21000

        Coordinates for packages such as GAMESS and MOPAC are also supported:

        >>> water_dimer = atoms.read_pybel("data/water-dimer.xyz")
        >>> print(water_dimer.to_string("gamin"))
        O      8.0     -1.6289300000   -0.0413800000    0.3713700000
        H      1.0     -0.6980300000   -0.0916800000    0.0933700000
        H      1.0     -2.0666300000   -0.7349800000   -0.1366300000
        O      8.0      1.2145700000    0.0317200000   -0.2762300000
        H      1.0      1.4492700000    0.9167200000   -0.5857300000
        H      1.0      1.7297700000   -0.0803800000    0.5338700000
        >>> print(water_dimer.to_string("mop"))
        O  -1.62893 1 -0.04138 1  0.37137 1
        H  -0.69803 1 -0.09168 1  0.09337 1
        H  -2.06663 1 -0.73498 1 -0.13663 1
        O   1.21457 1  0.03172 1 -0.27623 1
        H   1.44927 1  0.91672 1 -0.58573 1
        H   1.72977 1 -0.08038 1  0.53387 1

        Constraining of cartesian coordinates works with MOPAC format:

        >>> print(water_dimer.to_string("mop", constraints=(0, 3)))
        O  -1.62893 0 -0.04138 0  0.37137 0
        H  -0.69803 1 -0.09168 1  0.09337 1
        H  -2.06663 1 -0.73498 1 -0.13663 1
        O   1.21457 0  0.03172 0 -0.27623 0
        H   1.44927 1  0.91672 1 -0.58573 1
        H   1.72977 1 -0.08038 1  0.53387 1

        Fragment identification is supported for ``"ADF"`` and ``"ORCA"``
        dialects:

        >>> print(water_dimer.to_string("xyz", dialect="ADF",
        ...                             fragment_id="dimer"))
        O         -1.62893       -0.04138        0.37137       f=dimer
        H         -0.69803       -0.09168        0.09337       f=dimer
        H         -2.06663       -0.73498       -0.13663       f=dimer
        O          1.21457        0.03172       -0.27623       f=dimer
        H          1.44927        0.91672       -0.58573       f=dimer
        H          1.72977       -0.08038        0.53387       f=dimer
        >>> print(water_dimer.to_string("xyz", dialect="ORCA",
        ...                             fragment_id=1))
        O(1)      -1.62893       -0.04138        0.37137
        H(1)      -0.69803       -0.09168        0.09337
        H(1)      -2.06663       -0.73498       -0.13663
        O(1)       1.21457        0.03172       -0.27623
        H(1)       1.44927        0.91672       -0.58573
        H(1)       1.72977       -0.08038        0.53387

        """
        s = self.to_pybel().write(format).strip()

        if dialect is None:
            dialect = "standard"
        dialect = dialect.lower()

        if format == "xyz":
            natom, comment, body = s.split("\n", 2)

            if dialect in {"adf", "orca", "standard"}:
                if fragment_id is not None:
                    if dialect == "adf":
                        body = \
                            "\n".join(["{}       f={}".format(line,
                                                              fragment_id)
                                       for line in body.split("\n")])
                    elif dialect == "orca":
                        fragment_id = "({})".format(fragment_id)
                        body = \
                            "\n".join([line.replace(" " * len(fragment_id),
                                                    fragment_id, 1)
                                       for line in body.split("\n")])
                    else:
                        raise KeyError
            else:
                raise KeyError

            if with_header:
                s = "\n".join([natom, comment, body])
            else:
                s = body

        elif format == "gamin":
            lines = s.split("\n")
            begin = "\n".join([line.strip() for line in lines[:5]])
            body = "\n".join([line.strip() for line in lines[5:-1]])

            if with_header:
                s = "\n".join([begin, body])
            else:
                s = body

        elif format == "mop":
            chunks = s.split("\n", 2)
            begin = "\n".join([line.strip() for line in chunks[:2]])
            body = chunks[2].strip()

            if constraints is not None:
                body = body.split("\n")
                for i in constraints:
                    body[i] = _re.sub(' 1( |$)', ' 0\g<1>', body[i])
                body = "\n".join(body)

            if with_header:
                s = "\n".join([begin, body])
            else:
                s = body

        return s.strip()


def read_cclib(value, name=None):
    """
    Create an `Atoms` object from data attributes parsed by cclib.

    `cclib <https://cclib.github.io/>`_ is an open source library, written in
    Python, for parsing and interpreting the results (logfiles) of
    computational chemistry packages.

    Parameters
    ----------
    value : `str`, `cclib.parser.logfileparser.Logfile`, `cclib.parser.data.ccData`
        A path to a logfile, or either a cclib job object (i.e., from
        `cclib.ccopen`), or cclib data object (i.e., from ``job.parse()``).
    name : `str`, optional
        Name for chemical species. If not given, this is set to the logfile
        path, if known. Chemical equations mention this name when refering to
        the returned object.

    Returns
    -------
    molecule : `Atoms`
        All attributes obtainable by cclib are made available as attributes in
        the returned object.

    Examples
    --------
    >>> from pyrrole.atoms import read_cclib
    >>> molecule = read_cclib('data/pyrrole.out')
    >>> molecule.atomnos
    array([6, 6, 6, 6, 7, 1, 1, 1, 1, 1], dtype=int32)
    >>> molecule.charge
    0

    """
    if isinstance(value, _logfileparser.Logfile):
        # TODO: test this case.
        jobfilename = value.filename
        ccdata = value.parse()
    elif isinstance(value, _data.ccData):
        # TODO: test this case.
        jobfilename = None
        ccdata = value
    else:
        # TODO: test this case.
        ccobj = _cclib.ccopen(value)
        jobfilename = ccobj.filename
        ccdata = ccobj.parse()

    if name is None:
        name = jobfilename

    attributes = ccdata.getattributes()
    attributes.update({
        'name': name,
        'jobfilename': jobfilename,
    })

    return Atoms(attributes)


def read_pybel(value, name=None):
    """
    Create an `Atoms` object from content parsed by Pybel.

    `Pybel <https://openbabel.org/docs/dev/UseTheLibrary/Python_Pybel.html>`_
    is a Python module that simplifies access to the OpenBabel API, a chemical
    toolbox designed to speak the many languages of chemical data. It’s an
    open, collaborative project allowing anyone to search, convert, analyze, or
    store data from molecular modeling, chemistry, solid-state materials,
    biochemistry, and related areas.

    Parameters
    ----------
    value : `str`, `pybel.Molecule`, `openbabel.OBMol`
        A path to a file, or either a Pybel Molecule object, or OpenBabel
        OBMol.
    name : `str`, optional
        Name for chemical species. If not given, this is set to the file path,
        if known. Chemical equations mention this name when refering to the
        returned object.

    Returns
    -------
    molecule : `Atoms`
        All attributes convertible from Pybel to cclib are made available as
        attributes in the returned object.

    Notes
    -----
    The following attributes are converted from Pybel to cclib: `atomcoords`,
    `atommasses`, `atomnos`, `natom`, `charge` and `mult`. One must keep in
    mind that `charge` and `mult` are not always reliable, since these are
    often calculated from atomic formal charges.

    Examples
    --------
    >>> from pyrrole.atoms import read_pybel
    >>> molecule = read_pybel('data/pyrrole.xyz')
    >>> molecule.atomnos
    array([6, 6, 6, 6, 7, 1, 1, 1, 1, 1], dtype=int32)
    >>> molecule.natom
    10
    >>> molecule.charge
    0

    """
    if isinstance(value, _pb.Molecule):
        # TODO: test this case.
        jobfilename = None
        charge, mult = value.charge, value.spin
        ccdata = _makecclib(value.OBMol)
    elif isinstance(value, _ob.OBMol):
        # TODO: test this case.
        jobfilename = None
        charge, mult = value.GetTotalCharge(), value.GetTotalSpinMultiplicity()
        ccdata = _makecclib(value)
    else:
        # TODO: test this case.
        jobfilename = value
        _, jobfilename_ext = _os.path.splitext(jobfilename)

        # TODO: This only reads first structure.
        mol = next(_pb.readfile(jobfilename_ext[1:], jobfilename))
        charge, mult = mol.charge, mol.spin
        ccdata = _makecclib(mol.OBMol)

    if name is None:
        name = jobfilename

    attributes = ccdata.getattributes()
    attributes.update({
        'name': name,
        'jobfilename': jobfilename,
        'charge': charge,
        'mult': mult
    })

    return Atoms(attributes)


def create_data(*args):
    """
    Produce a single data object from an arbitrary number of different objects.

    This function returns a single `pandas.DataFrame` object from a collection
    of `Atoms` and `pandas.DataFrame` objects. The returned object, already
    indexed by `Atoms.name`, can be promptly used by e.g. `ChemicalSystem`.

    Parameters
    ----------
    *args : `pandas.DataFrame` or `Atoms`-like
        All positional arguments are assumed to be sources of data.

        `Atoms`-like objects (i.e. any object accepted by the `Atoms`
        constructor) become single row records in the final returned
        data object. `pandas.DataFrame` data table objects, on the other hand,
        are concatenated together (by using `pandas.DataFrame.concat`).

    Returns
    -------
    dataframe : `pandas.DataFrame`
        Resulting tabular data object. The returned object is guaranteed to be
        indexed by `Atoms.name`; if no column with this name exists at
        indexing time, a new column (with `None` values) is created for the
        purpose of indexing.

    Notes
    -----
    The returned `pandas.DataFrame` will be indexed by `Atoms.name` (see
    examples below), which might be the same as `Atoms.jobfilename` if no name
    was given to the constructor of `Atoms` (e.g. mapping).

    Examples
    --------
    >>> from pyrrole.atoms import Atoms, create_data, read_cclib
    >>> pyrrole = read_cclib('data/pyrrolate/pyrrole.out', 'pyrrole')
    >>> pyrrolate = read_cclib('data/pyrrolate/pyrrolate.out')
    >>> data = create_data(pyrrole, pyrrolate)
    >>> data['charge']
    name
    pyrrole                         0
    data/pyrrolate/pyrrolate.out   -1
    Name: charge, dtype: int64

    """
    def _prepare_data(data):
        if not isinstance(data, _pd.DataFrame):
            try:
                data = _pd.DataFrame([data.to_series()])
            except AttributeError:
                data = _pd.DataFrame([Atoms(data).to_series()])
        if data.index.name != "name":
            if "name" not in data.columns:
                data["name"] = None
            data = data.set_index("name")
        return data.reset_index()
    args = map(_prepare_data, args)

    dataframe = _pd.concat(args, sort=False)
    return dataframe.set_index("name")
