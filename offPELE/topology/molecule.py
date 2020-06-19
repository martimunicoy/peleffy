
# Global imports

# Local imports
from offPELE.template import Impact


class Atom(object):
    def __init__(self, _id=-1, OPLS_type=None, PDB_name=None, unknown=None,
                 z_matrix_x=None, z_matrix_y=None, z_matrix_z=None,
                 sigma=None, epsilon=None, charge=None, born_radius=None,
                 SASA_radius=None, nonpolar_gamma=None, nonpolar_alpha=None):
        self.id = _id
        self.OPLS_type = OPLS_type
        self.PDB_name = PDB_name
        self.unknown = unknown
        self.z_matrix_x = z_matrix_x
        self.z_matrix_y = z_matrix_y
        self.z_matrix_z = z_matrix_z
        self.sigma = sigma
        self.epsilon = epsilon
        self.charge = charge
        self.born_radius = born_radius  # Rad. Non Polar SGB
        self.SASA_radius = SASA_radius  # Rad. Non Polar Type
        self.nonpolar_gamma = nonpolar_gamma  # SGB Non Polar gamma
        self.nonpolar_alpha = nonpolar_alpha  # SGB Non Polar type


class Molecule(object):
    def __init__(self, path=None):
        if isinstance(path, str):
            from pathlib import Path
            extension = Path(path).suffix
            extension = extension.strip('.')
            if extension == 'pdb':
                self._initialize_from_pdb(path)
            else:
                raise ValueError(
                    '{} is not a valid extension'.format(extension))
        else:
            self._initialize()

    def _initialize(self):
        self._name = ''
        self._forcefield = ''
        self._atoms = list()
        self._bonds = list()
        self._angles = list()
        self._dihedrals = list()
        self._impropers = list()
        self._rdkit_molecule = None
        self._off_molecule = None

    def _initialize_from_pdb(self, path):
        self._initialize()
        print(' - Loading molecule from RDKit')

        try:
            from rdkit import Chem
        except ImportError:
            raise Exception('RDKit Python API not found')

        self._rdkit_molecule = Chem.rdmolfiles.MolFromPDBFile(path,
                                                              removeHs=False)

        try:
            from openforcefield.topology.molecule import Molecule as \
                OFFMolecule
        except ImportError:
            raise Exception('OpenForceField package not found')

        self._off_molecule = OFFMolecule.from_rdkit(self.rdkit_molecule)

    def parameterize(self, forcefield):
        try:
            from openforcefield.topology.molecule import Molecule as \
                OFFMolecule
            from openforcefield.typing.engines.smirnoff import ForceField as \
                OFFForceField
        except ImportError:
            raise Exception('OpenForceField package not found')

        if not isinstance(self._off_molecule, OFFMolecule):
            raise Exception('OpenForceField molecule was not initialized '
                            + 'correctly')

        print(' - Loading forcefield')
        if isinstance(forcefield, str):
            forcefield = OFFForceField(forcefield)
        elif isinstance(forcefield, OFFForceField):
            pass
        else:
            raise Exception('Invalid forcefield type')

        print(' - Computing partial charges with am1bcc')
        self._off_molecule.compute_partial_charges_am1bcc()

        self._build_atoms()

    def _assert_parameterized(self):
        try:
            assert self.off_molecule is not None
        except AssertionError:
            raise Exception('Molecule not parameterized')

    def _build_atoms(self):
        self._assert_parameterized()

        for atom in self.off_molecule.atoms:
            self._add_atom(Atom())

    def _add_atom(self, atom):
        self._atoms.append(atom)

    def get_pdb_atom_names(self):
        self._assert_parameterized()

        pdb_atom_names = list()

        for atom in self.rdkit_molecule.GetAtoms():
            pdb_info = atom.GetPDBResidueInfo()
            pdb_atom_names.append(pdb_info.GetName())

        return pdb_atom_names

    def to_impact(self, path):
        pass

    @property
    def off_molecule(self):
        return self._off_molecule

    @property
    def rdkit_molecule(self):
        return self._rdkit_molecule
