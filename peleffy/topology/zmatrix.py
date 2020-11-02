"""
This module contains the classes and functions to generate the zmatrix of
a molecule.
"""

from copy import deepcopy

import numpy as np

from peleffy.topology.molecule import DummyAtom


class ZMatrix(np.ndarray):
    """
    It generates the zmatrix of a molecule as a numpy.array.

    Inspired by the PlopRotTemp algorithm.
    """

    def __init__(self, molecule):
        """
        It initializes a ZMatrix object.

        Parameters
        ----------
        molecule : An peleffy.topology.Molecule
            A Molecule object to be written as an Impact file
        """
        # We will work on a copy of the molecule's object to modify it freely
        self._molecule = deepcopy(molecule)

    def __new__(cls, molecule):
        """
        It customizes the creation of the numpy.array.

        Parameters
        ----------
        molecule : An peleffy.topology.Molecule
            A Molecule object to be written as an Impact file
        """
        molecule = deepcopy(molecule)
        obj = np.zeros((len(molecule.atoms), 3)).view(cls)
        coords = cls._extract_coords(molecule)
        obj = cls._build_zmatrix(cls, obj, coords, molecule)

        return obj

    @staticmethod
    def _extract_coords(molecule):
        """
        It extracts the coordinates of the molecule's atoms.

        Parameters
        ----------
        molecule : An peleffy.topology.Molecule
            A Molecule object to be written as an Impact file

        Returns
        -------
        coords : list[tuple[float]]
            The coordinates of the molecule
        """
        coords = list()
        for atom in molecule.atoms:
            coords.append((atom.x, atom.y, atom.z))

        return coords

    @staticmethod
    def _get_absolute_parent(molecule):
        """
        It returns the absolute parent in the topology of the molecule.

        Parameters
        ----------
        molecule : An peleffy.topology.Molecule
            A Molecule object to be written as an Impact file

        Returns
        -------
        absolute_parent : peleffy.topology.molecule.Atom
            The absolute parent of the molecule
        """
        absolute_parent = list()
        for atom in molecule.atoms:
            if atom.parent is None:
                absolute_parent.append(atom)

        assert len(absolute_parent) == 1, 'Only 1 absolute parent is expected'

        return absolute_parent[0]

    @staticmethod
    def _calculate_bond(x1, y1, z1, x2, y2, z2):
        """
        It calculates the bond distance between two sets of coordinates.

        Parameters
        ----------
        x1 : float
            The x1 coordinate
        y1 : float
            The y1 coordinate
        z1 : float
            The z1 coordinate
        x2 : float
            The x2 coordinate
        y2 : float
            The y2 coordinate
        z2 : float
            The z2 coordinate

        Returns
        -------
        distance : float
            The bond distance between the two sets of coordinates
        """
        dx = x1 - x2
        dy = y1 - y2
        dz = z1 - z2

        return np.sqrt(dx * dx + dy * dy + dz * dz)

    @staticmethod
    def _calculate_angle(x1, y1, z1, x2, y2, z2, x3, y3, z3):
        """
        It calculates the angle between three sets of coordinates.

        Parameters
        ----------
        x1 : float
            The x1 coordinate
        y1 : float
            The y1 coordinate
        z1 : float
            The z1 coordinate
        x2 : float
            The x2 coordinate
        y2 : float
            The y2 coordinate
        z2 : float
            The z2 coordinate
        x3 : float
            The x3 coordinate
        y3 : float
            The y3 coordinate
        z3 : float
            The z3 coordinate

        Returns
        -------
        angle : float
            The angle between the three sets of coordinates
        """
        dx_12 = x1 - x2
        dy_12 = y1 - y2
        dz_12 = z1 - z2

        dx_31 = x3 - x1
        dy_31 = y3 - y1
        dz_31 = z3 - z1

        vdot = dx_12 * dx_31 + dy_12 * dy_31 + dz_12 * dz_31

        assert np.fabs(vdot) > 0, 'A non-zero angle is expected'

        d12 = np.sqrt(dx_12 * dx_12 + dy_12 * dy_12 + dz_12 * dz_12)
        d31 = np.sqrt(dx_31 * dx_31 + dy_31 * dy_31 + dz_31 * dz_31)

        xang = vdot / (d12 * d31)

        if xang - 1.0 > -0.0000000001:
            return 0.0
        elif xang + 1.0 < 0.0000000001:
            return np.pi

        return np.arccos(xang) * 180.0 / np.pi

    @staticmethod
    def _calculate_dihedral(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4):
        """
        It calculates the dihedral between four sets of coordinates.

        Parameters
        ----------
        x1 : float
            The x1 coordinate
        y1 : float
            The y1 coordinate
        z1 : float
            The z1 coordinate
        x2 : float
            The x2 coordinate
        y2 : float
            The y2 coordinate
        z2 : float
            The z2 coordinate
        x3 : float
            The x3 coordinate
        y3 : float
            The y3 coordinate
        z3 : float
            The z3 coordinate
        x4 : float
            The x4 coordinate
        y4 : float
            The y4 coordinate
        z4 : float
            The z4 coordinate

        Returns
        -------
        dihedral : float
            The dihedral between the four sets of coordinates
        """
        dx_12 = x1 - x2
        dy_12 = y1 - y2
        dz_12 = z1 - z2

        dx_32 = x3 - x2
        dy_32 = y3 - y2
        dz_32 = z3 - z2

        dx_34 = x3 - x4
        dy_34 = y3 - y4
        dz_34 = z3 - z4

        ax = dy_12 * dz_32 - dz_12 * dy_32
        ay = dz_12 * dx_32 - dx_12 * dz_32
        az = dx_12 * dy_32 - dy_12 * dx_32
        cx = dy_32 * dz_34 - dz_32 * dy_34
        cy = dz_32 * dx_34 - dx_32 * dz_34
        cz = dx_32 * dy_34 - dy_32 * dx_34

        rac = ax * cx + ay * cy + az * cz
        ra = ax * ax + ay * ay + az * az
        rc = cx * cx + cy * cy + cz * cz

        cosang = rac / np.sqrt(ra * rc)

        if cosang - 1.0 > -0.00000000001:
            phi = 0.0
        elif cosang + 1.0 < 0.00000000001:
            phi = np.pi
        else:
            phi = np.arccos(cosang)

        s = dx_12 * cx + dy_12 * cy + dz_12 * cz
        if (s < 0):
            phi = -phi  # to account for phi between pi and 2pi

        return phi * 180.0 / np.pi

    @staticmethod
    def _build_zmatrix(cls, obj, coords, molecule):
        """
        It buils the zmatrix.

        Parameters
        ----------
        cls : ZMatrix class
            The ZMatrix class
        obj : ZMatrix object
            The ZMatrix object
        coords : list[tuple[float]]
            The coordinates of the molecule
        molecule : An peleffy.topology.Molecule
            A Molecule object to be written as an Impact file

        Returns
        -------
        obj :  ZMatrix object
            The ZMatrix object which is a numpy.array with the corresponding
            zmatrix initialized.
        """
        dummy1 = DummyAtom(index=-3, PDB_name='DUM1', parent=None)
        dummy2 = DummyAtom(index=-2, PDB_name='DUM2', parent=dummy1)
        dummy3 = DummyAtom(index=-1, PDB_name='DUM3', parent=dummy2)
        absolute_parent = cls._get_absolute_parent(molecule)
        absolute_parent.set_parent(dummy3)

        dummy1.set_coords([0.0, 0.0, 0.0])  # at origin
        dummy2.set_coords([0.0, 0.0, 1.0])  # at z axis
        dummy3.set_coords([1.0, 0.0, 0.0])  # at x axis

        for i, atom in enumerate(molecule.atoms):
            atom1 = atom
            atom2 = atom1.parent
            atom3 = atom2.parent
            atom4 = atom3.parent

            assert atom1 is not None and atom2 is not None \
                and atom3 is not None and atom4 is not None, \
                'A None as parent is not expected'

            x1, y1, z1 = (atom1.x, atom1.y, atom1.z)
            x2, y2, z2 = (atom2.x, atom2.y, atom2.z)
            x3, y3, z3 = (atom3.x, atom3.y, atom3.z)
            x4, y4, z4 = (atom4.x, atom4.y, atom4.z)

            obj[i][0] = cls._calculate_bond(x1, y1, z1,
                                            x2, y2, z2)
            obj[i][1] = cls._calculate_angle(x1, y1, z1,
                                             x2, y2, z2,
                                             x3, y3, z3)
            obj[i][2] = cls._calculate_dihedral(x1, y1, z1,
                                                x2, y2, z2,
                                                x3, y3, z3,
                                                x4, y4, z4)

        return obj
