"""
This module defines the basic topology elements of a molecule.
"""


from simtk import unit


class _TopologyElement(object):
    """
    A wrapper for any topological element.
    """

    _name = None
    _writable_attrs = []

    class TopologyIterator(object):
        """
        An iterator for topological elements that iterates over their
        attributes in an ordered way.

        It is useful when writing topological elements to file.
        """

        def __init__(self, top_el):
            """
            It initiates a TopologyIterator object.

            Parameters
            ----------
            top_el : a TopologyElement object
                The topology element to iterate on.
            """
            self._index = int(0)
            self._top_el = top_el

        def __next__(self):
            """
            It returns the next item for the iteration.

            Returns
            -------
            attr_name : str
                The name of the attribute
            attr_value : float
                The value of the attribute
            """
            if self._index == len(self._top_el._writable_attrs):
                raise StopIteration

            attr_name = self._top_el._writable_attrs[self._index]
            attr_value = getattr(self._top_el, attr_name)
            self._index += 1

            return attr_name, attr_value

    @property
    def name(self):
        """
        The name that this topological element has.

        Returns
        -------
        name : str
            The name of the topological element
        """
        return self._name

    @property
    def n_writable_attrs(self):
        """
        The number of writable attributes this topological element has.

        Returns
        -------
        n_writable_attrs : int
            The number of writable attributes
        """
        return len(self._writable_attrs)

    def __iter__(self):
        """
        It returns an instance of the TopologyIterator.

        Returns
        -------
        iterator : a TopologyIterator
            The TopologyIterator object
        """
        return self.TopologyIterator(self)

    def __repr__(self):
        """
        It returns the representation string of this topological element.

        Returns
        -------
        repr_string : str
            The representation string
        """
        repr_string = '{}('.format(self._name)
        attrs = [attr for attr in self]
        for attr_name, value in attrs[:-1]:
            repr_string += '{}={}, '.format(attr_name, value)
        repr_string += '{}={})'.format(*attrs[-1])

        return repr_string

    def __str__(self):
        """
        It returns the readable representation string of this topological
        element.

        Returns
        -------
        str_string : str
            The readable representation string
        """
        return self.__repr__()


class Atom(_TopologyElement):
    """
    It represents the properties of an atom.
    """

    _name = 'Atom'
    _writable_attrs = ['index', 'PDB_name', 'OPLS_type']

    def __init__(self, index=-1, core=None, OPLS_type=None, PDB_name=None,
                 unknown=None, x=None, y=None, z=None, sigma=None,
                 epsilon=None, charge=None, born_radius=None, SASA_radius=None,
                 nonpolar_gamma=None, nonpolar_alpha=None, parent=None):
        """
        It initializes an Atom object.

        Parameters
        ----------
        index : int
            The index of the atom
        core : bool
            Whether atom is in the core or in a branch
        OPLS_type : str
            The OPLS type of the atom
        PDB_name : str
            The PDB name of the atom
        unknown : int
            The unknown value of the atom
        x : float
            The x coordinate of the atom
        y : float
            The y coordinate of the atom
        z : float
            The z coordinate of the atom
        sigma : simtk.unit.Quantity
            The sigma parameter of the atom
        epsilon : simtk.unit.Quantity
            The epsilon parameter of the atom
        charge : simtk.unit.Quantity
            The partial charge parameter of the atom
        born_radius : simtk.unit.Quantity
            The SGB born parameter radius of the atom
        SASA_radius : simtk.unit.Quantity
            The SASA radius parameter of the atom
        nonpolar_gamma : simtk.unit.Quantity
            The nonpolar gamma parameter of the atom
        nonpolar_alpha : simtk.unit.Quantity
            The nonpolar alpha parameter of the atom
        parent : peleffy.topology.Atom
            The parent of the atom
        """
        self._index = index
        self._core = core
        self._OPLS_type = OPLS_type
        self._PDB_name = PDB_name
        self._unknown = unknown
        self._x = x
        self._y = y
        self._z = z
        self._sigma = sigma
        self._epsilon = epsilon
        self._charge = charge
        self._born_radius = born_radius  # Rad. Non Polar SGB
        self._SASA_radius = SASA_radius  # Rad. Non Polar Type
        self._nonpolar_gamma = nonpolar_gamma  # SGB Non Polar gamma
        self._nonpolar_alpha = nonpolar_alpha  # SGB Non Polar type
        self._parent = parent

    def set_index(self, index):
        """
        It sets the index of the atom.

        Parameters
        ----------
        index : int
            The index of this Atom object
        """
        self._index = index

    def set_as_core(self):
        """It sets the atom as core"""
        self._core = True

    def set_as_branch(self):
        """It sets the atom as branch"""
        self._core = False

    def set_parent(self, parent):
        """
        It sets the parent of the atom.

        Parameters
        ----------
        parent : peleffy.topology.Atom
            The parent of the atom
        """
        self._parent = parent

    def set_coords(self, coords):
        """
        It sets the coordinates of the atom.

        Parameters
        ----------
        coords : list
            The coordinates array to set to this Atom object
        """
        assert len(coords) == 3, '3D array is expected'

        self._x, self._y, self._z = coords

    def set_OPLS_type(self, OPLS_type):
        """
        It sets the OPLS type of the atom.

        Parameters
        ----------
        OPLS_type : str
            The OPLS type to set to this Atom object
        """
        self._OPLS_type = OPLS_type

    def set_sigma(self, sigma):
        """
        It sets the sigma of the atom.

        Parameters
        ----------
        sigma : float
            The sigma to set to this Atom object
        """
        self._sigma = sigma

    def set_epsilon(self, epsilon):
        """
        It sets the epsilon of the atom.

        Parameters
        ----------
        epsilon : float
            The epsilon to set to this Atom object
        """
        self._epsilon = epsilon

    def set_charge(self, charge):
        """
        It sets the charge of the atom.

        Parameters
        ----------
        charge : float
            The charge to set to this Atom object
        """
        self._charge = charge

    def set_born_radius(self, born_radius):
        """
        It sets the Born radius of the atom.

        Parameters
        ----------
        born_radius : float
            The Born radius to set to this Atom object
        """
        self._born_radius = born_radius

    def set_SASA_radius(self, SASA_radius):
        """
        It sets the SASA radius of the atom.

        Parameters
        ----------
        SASA_radius : float
            The SASA radius to set to this Atom object
        """
        self._SASA_radius = SASA_radius

    def set_nonpolar_gamma(self, nonpolar_gamma):
        """
        It sets the nonpolar gamma of the atom.

        Parameters
        ----------
        nonpolar_gamma : float
            The nonpolar gamma to set to this Atom object
        """
        self._nonpolar_gamma = nonpolar_gamma

    def set_nonpolar_alpha(self, nonpolar_alpha):
        """
        It sets the nonpolar alpha of the atom.

        Parameters
        ----------
        nonpolar_alpha : float
            The nonpolar alpha to set to this Atom object
        """
        self._nonpolar_alpha = nonpolar_alpha

    @property
    def index(self):
        """
        Atom's index.

        Returns
        -------
        index : int
            The index of this Atom object
        """
        return self._index

    @property
    def core(self):
        """
        Atom's core position.

        Returns
        -------
        core : bool
            Whether this Atom object is in the core or in a branch
        """
        return self._core

    @property
    def OPLS_type(self):
        """
        Atom's OPLS type.

        .. todo ::

           * Consider removing any reference to OPLS, if possible
             Otherwise, use SMIRks to find the best match

        Returns
        -------
        OPLS_type : str
            The OLPS type of this Atom object
        """
        return self._OPLS_type

    @property
    def PDB_name(self):
        """
        Atom's PDB name.

        .. todo ::

           * Consider removing any reference to OPLS, if possible
             Otherwise, use SMIRks to find the best match

        Returns
        -------
        PDB_name : str
            The PDB name of this Atom object
        """
        return self._PDB_name

    @property
    def unknown(self):
        """
        Atom's unknown int.

        .. todo ::

           * Review the actual purpose of this attribute in PELE

        Returns
        -------
        unknown : int
            The unknown int of this Atom object
        """
        return self._unknown

    @property
    def x(self):
        """
        Atom's x coordinate.

        Returns
        -------
        x : float
            The x coordinate of this Atom object
        """
        return self._x

    @property
    def y(self):
        """
        Atom's y coordinate.

        Returns
        -------
        y : float
            The y coordinate of this Atom object
        """
        return self._y

    @property
    def z(self):
        """
        Atom's z coordinate.

        Returns
        -------
        z : float
            The z coordinate of this Atom object
        """
        return self._z

    @property
    def sigma(self):
        """
        Atom's sigma.

        Returns
        -------
        sigma : simtk.unit.Quantity
            The sigma parameter of this Atom object
        """
        return self._sigma

    @property
    def epsilon(self):
        """
        Atom's epsilon.

        Returns
        -------
        epsilon : simtk.unit.Quantity
            The epsilon parameter of this Atom object
        """
        return self._epsilon

    @property
    def charge(self):
        """
        Atom's charge.

        Returns
        -------
        charge : simtk.unit.Quantity
            The charge parameter of this Atom object
        """
        return self._charge

    @property
    def born_radius(self):
        """
        Atom's born radius.

        Returns
        -------
        born_radius : simtk.unit.Quantity
            The SGB Born radius parameter of this Atom object
        """
        return self._born_radius

    @property
    def SASA_radius(self):
        """
        Atom's SASA radius.

        Returns
        -------
        SASA_radius : simtk.unit.Quantity
            The SASA radius parameter of this Atom object
        """
        return self._SASA_radius

    @property
    def nonpolar_gamma(self):
        """
        Atom's nonpolar gamma.

        Returns
        -------
        nonpolar_gamma : simtk.unit.Quantity
            The nonpolar gamma parameter of this Atom object
        """
        return self._nonpolar_gamma

    @property
    def nonpolar_alpha(self):
        """
        Atom's nonpolar alpha.

        Returns
        -------
        nonpolar_alpha : simtk.unit.Quantity
            The nonpolar alpha parameter of this Atom object
        """
        return self._nonpolar_alpha

    @property
    def parent(self):
        """
        Atom's parent.

        Returns
        -------
        parent : simtk.unit.Quantity
            The nonpolar gamma parameter of this Atom object
        """
        return self._parent


class DummyAtom(Atom):
    """
    It represents a dummy atom.
    """

    def __init__(self, index=-1, PDB_name='DUMM', parent=None):
        """
        It initializes a DummyAtom object.

        Parameters
        ----------
        index : int
            The index of the atom
        PDB_name : str
            The PDB name of the atom
        parent : peleffy.topology.Atom
            The parent of the atom
        """
        if parent is None:
            parent = self
        super().__init__(index, False, None, PDB_name, None, None, None, None,
                         None, None, None, None, None, None, None, parent)


class Bond(_TopologyElement):
    """
    It represents a bond in the topology.
    """

    _name = 'Bond'
    _writable_attrs = ['atom1_idx', 'atom2_idx', 'spring_constant', 'eq_dist']

    def __init__(self, index=-1, atom1_idx=None, atom2_idx=None,
                 spring_constant=None, eq_dist=None):
        """
        It initiates a Bond object.

        Parameters
        ----------
        index : int
            The index of this Bond object
        atom1_idx : int
            The index of the first atom involved in this Bond
        atom2_idx : int
            The index of the second atom involved in this Bond
        spring_constant : simtk.unit.Quantity
            The spring constant of this Bond
        eq_dist : simtk.unit.Quantity
            The equilibrium distance of this Bond
        """
        self._index = index
        self._atom1_idx = atom1_idx
        self._atom2_idx = atom2_idx
        self._spring_constant = spring_constant
        self._eq_dist = eq_dist

    def set_atom1_idx(self, index):
        """
        It sets atom1's index.

        Parameters
        ----------
        index : int
            The index of the first atom involved in this Bond
        """
        self._atom1_idx = index

    def set_atom2_idx(self, index):
        """
        It sets atom2's index.

        Parameters
        ----------
        index : int
            The index of the second atom involved in this Bond
        """
        self._atom2_idx = index

    @property
    def index(self):
        """
        Bond's index.

        Returns
        -------
        index : int
            The index of this Bond object
        """
        return self._index

    @property
    def atom1_idx(self):
        """
        Bond's atom1 index.

        Returns
        -------
        atom1_idx : int
            The index of the first atom involved in this Bond object
        """
        return self._atom1_idx

    @property
    def atom2_idx(self):
        """
        Bond's atom2 index.

        Returns
        -------
        atom2_idx : int
            The index of the second atom involved in this Bond object
        """
        return self._atom2_idx

    @property
    def spring_constant(self):
        """
        Bond's spring constant.

        Returns
        -------
        spring_constant : simtk.unit.Quantity
            The spring constant of this Bond object
        """
        return self._spring_constant

    @property
    def eq_dist(self):
        """
        Bond's equilibrium distance.

        Returns
        -------
        eq_dist : simtk.unit.Quantity
            The equilibrium distance of this Bond object
        """
        return self._eq_dist


class Angle(_TopologyElement):
    """
    It represents an angle in the topology.
    """

    _name = 'Angle'
    _writable_attrs = ['atom1_idx', 'atom2_idx', 'atom3_idx',
                       'spring_constant', 'eq_angle']

    def __init__(self, index=-1, atom1_idx=None, atom2_idx=None,
                 atom3_idx=None, spring_constant=None, eq_angle=None):
        """
        It initiates an Angle object.

        Parameters
        ----------
        index : int
            The index of this Angle object
        atom1_idx : int
            The index of the first atom involved in this Angle
        atom2_idx : int
            The index of the second atom involved in this Angle
        atom3_idx : int
            The index of the third atom involved in this Angle
        spring_constant : simtk.unit.Quantity
            The spring constant of this Angle
        eq_angle : simtk.unit.Quantity
            The equilibrium angle of this Angle
        """
        self._index = index
        self._atom1_idx = atom1_idx
        self._atom2_idx = atom2_idx
        self._atom3_idx = atom3_idx
        self._spring_constant = spring_constant
        self._eq_angle = eq_angle

    def set_atom1_idx(self, index):
        """
        It sets atom1's index.

        Parameters
        ----------
        index : int
            The index of the first atom involved in this Angle
        """
        self._atom1_idx = index

    def set_atom2_idx(self, index):
        """
        It sets atom2's index.

        Parameters
        ----------
        index : int
            The index of the second atom involved in this Angle
        """
        self._atom2_idx = index

    def set_atom3_idx(self, index):
        """
        It sets atom3's index.

        Parameters
        ----------
        index : int
            The index of the third atom involved in this Angle
        """
        self._atom3_idx = index

    @property
    def index(self):
        """
        Angle's index.

        Returns
        -------
        index : int
            The index of this Angle object
        """
        return self._index

    @property
    def atom1_idx(self):
        """
        Angle's atom1 index.

        Returns
        -------
        atom1_idx : int
            The index of the first atom involved in this Angle object
        """
        return self._atom1_idx

    @property
    def atom2_idx(self):
        """
        Angle's atom2 index.

        Returns
        -------
        atom1_idx : int
            The index of the second atom involved in this Angle object
        """
        return self._atom2_idx

    @property
    def atom3_idx(self):
        """
        Angle's atom3 index.

        Returns
        -------
        atom1_idx : int
            The index of the third atom involved in this Angle object
        """
        return self._atom3_idx

    @property
    def spring_constant(self):
        """
        Angle's spring constant.

        Returns
        -------
        spring_constant : simtk.unit.Quantity
            The spring constant of this Angle object
        """
        return self._spring_constant

    @property
    def eq_angle(self):
        """
        Angle's equilibrium angle.

        Returns
        -------
        eq_angle : simtk.unit.Quantity
            The equilibrium angle of this Angle object
        """
        return self._eq_angle


class Dihedral(_TopologyElement):
    """
    It represents a dihedral in the topology.

    It can be a proper or an improper dihedral.
    """

    _name = 'Dihedral'

    def __init__(self, index=-1, atom1_idx=None, atom2_idx=None,
                 atom3_idx=None, atom4_idx=None, periodicity=None,
                 prefactor=None, constant=None, phase=None):
        """
        It initiates an Dihedral object.

        Parameters
        ----------
        index : int
            The index of this Dihedral object
        atom1_idx : int
            The index of the first atom involved in this Dihedral
        atom2_idx : int
            The index of the second atom involved in this Dihedral
        atom3_idx : int
            The index of the third atom involved in this Dihedral
        atom4_idx : int
            The index of the fourth atom involved in this Dihedral
        periodicity : int
            The periodicity of this Dihedral
        prefactor : int
            The prefactor of this Dihedral
        constant : simtk.unit.Quantity
            The constant of this Dihedral
        phase : simtk.unit.Quantity
            The phase constant of this Dihedral
        """
        self._index = index
        self._atom1_idx = atom1_idx
        self._atom2_idx = atom2_idx
        self._atom3_idx = atom3_idx
        self._atom4_idx = atom4_idx
        self._periodicity = periodicity
        self._prefactor = prefactor
        self._constant = constant
        self._phase = phase

    def set_atom1_idx(self, index):
        """
        It sets atom1's index.

        Parameters
        ----------
        index : int
            The index of the first atom involved in this Dihedral
        """
        self._atom1_idx = index

    def set_atom2_idx(self, index):
        """
        It sets atom2's index.

        Parameters
        ----------
        index : int
            The index of the second atom involved in this Dihedral
        """
        self._atom2_idx = index

    def set_atom3_idx(self, index):
        """
        It sets atom3's index.

        Parameters
        ----------
        index : int
            The index of the third atom involved in this Dihedral
        """
        self._atom3_idx = index

    def set_atom4_idx(self, index):
        """
        It sets atom4's index.

        Parameters
        ----------
        index : int
            The index of the fourth atom involved in this Dihedral
        """
        self._atom4_idx = index

    def plot(self):
        """
        It plots this Dihedral as a function of phi angle.
        """
        from matplotlib import pyplot
        import numpy as np

        x = unit.Quantity(np.arange(0, np.pi, 0.1), unit=unit.radians)
        pyplot.plot(x, self.constant * (1 + self.prefactor
                                        * np.cos(self.periodicity * x
                                                 + self.phase)),
                    'r--')

        pyplot.show()

    @property
    def index(self):
        """
        Dihedral's index.

        Returns
        -------
        index : int
            The index of this Dihedral object
        """
        return self._index

    @property
    def atom1_idx(self):
        """
        Dihedral's atom1 index.

        Returns
        -------
        atom1_idx : int
            The index of the first atom involved in this Dihedral object
        """
        return self._atom1_idx

    @property
    def atom2_idx(self):
        """
        Dihedral's atom2 index.

        Returns
        -------
        atom1_idx : int
            The index of the second atom involved in this Dihedral object
        """
        return self._atom2_idx

    @property
    def atom3_idx(self):
        """
        Dihedral's atom3 index.

        Returns
        -------
        atom1_idx : int
            The index of the third atom involved in this Dihedral object
        """
        return self._atom3_idx

    @property
    def atom4_idx(self):
        """
        Dihedral's atom4 index.

        Returns
        -------
        atom1_idx : int
            The index of the fourth atom involved in this Dihedral object
        """
        return self._atom4_idx

    @property
    def periodicity(self):
        """
        Dihedral's periodicity.

        Returns
        -------
        periodicity : int
            The periodicity this Dihedral object
        """
        return self._periodicity

    @property
    def prefactor(self):
        """
        Dihedral's prefactor.

        Returns
        -------
        prefactor : int
            The prefactor this Dihedral object
        """
        return self._prefactor

    @property
    def constant(self):
        """
        Dihedral's constant.

        Returns
        -------
        constant : unit.simtk.Quantity
            The constant of this Dihedral object
        """
        return self._constant

    @property
    def phase(self):
        """
        Dihedral's phase constant.

        Returns
        -------
        phase : unit.simtk.Quantity
            The phase constant of this Dihedral object
        """
        return self._phase


class Proper(Dihedral):
    """
    It represents a proper dihedral in the topology.
    """

    _name = 'Proper'
    exclude = False
    _writable_attrs = ['atom1_idx', 'atom2_idx', 'atom3_idx', 'atom4_idx',
                       'constant', 'prefactor', 'periodicity', 'phase']

    def include_in_14_list(self):
        """
        It includes this proper dihedral in PELE's 1-4 list.
        """
        self.exclude = False

    def exclude_from_14_list(self):
        """
        It excludes this proper dihedral from PELE's 1-4 list by
        setting the index of the third atom to negative.
        """
        self.exclude = True


class Improper(Dihedral):
    """
    It represents an improper dihedral in the topology.
    """

    _name = 'Improper'
    _writable_attrs = ['atom1_idx', 'atom2_idx', 'atom3_idx', 'atom4_idx',
                       'constant', 'prefactor', 'periodicity']


class OFFDihedral(_TopologyElement):
    """
    It represents a dihedral in the Open Force Field's topology.
    """

    _name = 'OFFDihedral'
    _writable_attrs = ['atom1_idx', 'atom2_idx', 'atom3_idx', 'atom4_idx',
                       'periodicity', 'phase', 'k', 'idivf']
    _to_PELE_class = Dihedral

    def __init__(self, index=-1, atom1_idx=None, atom2_idx=None,
                 atom3_idx=None, atom4_idx=None, periodicity=None,
                 phase=None, k=None, idivf=None):
        """
        It initiates an Dihedral object.

        Parameters
        ----------
        index : int
            The index of this Dihedral object
        atom1_idx : int
            The index of the first atom involved in this Dihedral
        atom2_idx : int
            The index of the second atom involved in this Dihedral
        atom3_idx : int
            The index of the third atom involved in this Dihedral
        atom4_idx : int
            The index of the fourth atom involved in this Dihedral
        periodicity : int
            The periodicity of this Dihedral
        phase : simtk.unit.Quantity
            The phase of this Dihedral
        k : simtk.unit.Quantity
            The constant of this Dihedral
        idivf : int
            The idivf of this Dihedral
        """
        self.index = index
        self.atom1_idx = atom1_idx
        self.atom2_idx = atom2_idx
        self.atom3_idx = atom3_idx
        self.atom4_idx = atom4_idx
        self.periodicity = periodicity
        self.phase = phase
        self.k = k
        self.idivf = idivf

    def _check_up(self):
        """
        It performs some parameter check ups.

        Raises
        ------
        AssertionError
            If any unexpected value is found
        """
        pass

    def to_PELE(self):
        """
        It converts this Open Force Field Dihedral object into a
        PELE-compatible one.

        .. todo ::

           * Fully cover all OpenFF dihedral parameters

        Returns
        -------
        PELE_dihedral : a Dihedral
            The PELE-compatible Dihedral object
        """
        if (self.periodicity is None or self.phase is None
                or self.k is None or self.idivf is None):
            return None

        try:
            self._check_up()
        except AssertionError as e:
            raise ValueError('Invalid value found: {}'.format(e))

        PELE_phase = self.phase

        if self.phase.value_in_unit(unit.degree) == 180:
            PELE_prefactor = -1
            PELE_phase = unit.Quantity(value=0.0, unit=unit.degree)
        else:
            PELE_prefactor = 1

        PELE_constant = self.k / self.idivf

        PELE_dihedral_kwargs = {'index': self.index,
                                'atom1_idx': self.atom1_idx,
                                'atom2_idx': self.atom2_idx,
                                'atom3_idx': self.atom3_idx,
                                'atom4_idx': self.atom4_idx,
                                'periodicity': self.periodicity,
                                'prefactor': PELE_prefactor,
                                'constant': PELE_constant,
                                'phase': PELE_phase}

        return self._to_PELE_class(**PELE_dihedral_kwargs)

    def plot(self):
        """
        It plots this Dihedral as a function of phi angle.
        """
        from matplotlib import pyplot
        import numpy as np

        x = unit.Quantity(np.arange(0, np.pi, 0.1), unit=unit.radians)
        pyplot.plot(x,
                    self.k * (1 + np.cos(self.periodicity * x - self.phase)),
                    'r--')

        pyplot.show()


class OFFProper(OFFDihedral):
    """
    It represents a proper dihedral in the Open Force Field's topology.
    """

    _name = 'OFFProper'
    _to_PELE_class = Proper

    def _check_up(self):
        """
        It performs some parameter check ups.

        .. todo ::

            * Periodicity can also equal 5 and currently it is not supported
             by PELE

        Raises
        ------
        AssertionError
            If any unexpected value is found
        """
        assert self.periodicity in (1, 2, 3, 4, 5, 6), 'Expected values ' \
            'for periodicity are 1, 2, 3, 4, 5 or 6, obtained ' \
            '{}'.format(self.periodicity)

        # proper's idivfs must always be 1 --> Apparently not: COC(O)OC, t84, t86
        # assert self.idivf == 1, 'The expected value for idivf is 1 ' \
        #     'for propers, obtained {}'.format(self.idivf)

        # The next version of PELE will be compatible with phase values
        # other than 0 and 180 degrees to fully cover all OpenFF dihedrals
        # assert self.phase.value_in_unit(unit.degree) in (0, 180), \
        #     'Expected values for phase are 0 or 180, obtained ' \
        #     '{}'.format(self.phase)


class OFFImproper(OFFDihedral):
    """
    It represents an improper dihedral in the Open Force Field's topology.
    """

    _name = 'OFFImproper'
    _to_PELE_class = Improper

    def _check_up(self):
        """
        It performs some parameter check ups.

        .. todo ::

            * Periodicity can also equal 5 and currently it is not supported
             by PELE

        Raises
        ------
        AssertionError
            If any unexpected value is found
        """
        assert self.periodicity in (1, 2, 3, 4, 5, 6), 'Expected values ' \
            'for periodicity are 1, 2, 3, 4, 5 or 6, obtained ' \
            '{}'.format(self.periodicity)

        assert self.phase.value_in_unit(unit.degree) in (0, 180), \
            'Expected values for phase are 0 or 180 in impropers, ' \
            'obtained {}'.format(self.phase)
