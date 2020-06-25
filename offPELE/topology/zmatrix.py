from copy import deepcopy

import numpy as np

from offPELE.topology.molecule import Atom, DummyAtom


class ZMatrix(np.ndarray):
    """
    Inspired by the PlopRotTemp algorithm
    """

    def __init__(self, molecule):
        # We will work on a copy of the molecule's object to modify it freely
        self._molecule = deepcopy(molecule)

    def __new__(cls, molecule):
        molecule = deepcopy(molecule)
        obj = np.zeros((len(molecule.atoms), 3)).view(cls)
        coords = cls._extract_coords(molecule)
        obj = cls._build_zmatrix(cls, obj, coords, molecule)

        return obj

    def _extract_coords(molecule):
        coords = list()
        for atom in molecule.atoms:
            coords.append((atom.x, atom.y, atom.z))

        return coords

    def _get_absolute_parent(molecule):
        absolute_parent = list()
        for atom in molecule.atoms:
            if atom.parent is None:
                absolute_parent.append(atom)

        assert len(absolute_parent) == 1, 'Only 1 absolute parent is expected'

        return absolute_parent[0]

    def _calculate_bond(x1, y1, z1, x2, y2, z2):
        dx = x1 - x2
        dy = y1 - y2
        dz = z1 - z2

        return np.sqrt(dx * dx + dy * dy + dz * dz)

    def _calculate_angle(x1, y1, z1, x2, y2, z2, x3, y3, z3):
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

    def _calculate_dihedral(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4):
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

    def _build_zmatrix(cls, obj, coords, molecule):

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
                'None parent is not expected'

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
