# Global imports
from simtk import unit


class TopologyElement(object):
    _name = None
    _writable_attrs = []

    class TopologyIterator(object):
        def __init__(self, top_el):
            self._index = int(0)
            self._top_el = top_el

        def __next__(self):
            if self._index == len(self._top_el._writable_attrs):
                raise StopIteration

            attr_name = self._top_el._writable_attrs[self._index]
            attr_value = getattr(self._top_el, attr_name)
            self._index += 1

            return attr_name, attr_value

    @property
    def name(self):
        return self._name

    @property
    def n_writable_attrs(self):
        return len(self._writable_attrs)

    def __iter__(self):
        return self.TopologyIterator(self)

    def __repr__(self):
        repr_string = '{}('.format(self._name)
        attrs = [attr for attr in self]
        for attr_name, value in attrs[:-1]:
            repr_string += '{}={}, '.format(attr_name, value)
        repr_string += '{}={})'.format(*attrs[-1])

        return repr_string

    def __str__(self):
        return self.__repr__()


class Bond(TopologyElement):
    _name = 'Bond'
    _writable_attrs = ['atom1_idx', 'atom2_idx', 'spring_constant', 'eq_dist']

    def __init__(self, index=-1, atom1_idx=None, atom2_idx=None,
                 spring_constant=None, eq_dist=None):
        self._index = index
        self._atom1_idx = atom1_idx
        self._atom2_idx = atom2_idx
        self._spring_constant = spring_constant
        self._eq_dist = eq_dist

    def set_atom1_idx(self, index):
        self._atom1_idx = index

    def set_atom2_idx(self, index):
        self._atom2_idx = index

    @property
    def index(self):
        return self._index

    @property
    def atom1_idx(self):
        return self._atom1_idx

    @property
    def atom2_idx(self):
        return self._atom2_idx

    @property
    def spring_constant(self):
        return self._spring_constant

    @property
    def eq_dist(self):
        return self._eq_dist


class Angle(TopologyElement):
    _name = 'Angle'
    _writable_attrs = ['atom1_idx', 'atom2_idx', 'atom3_idx',
                       'spring_constant', 'eq_angle']

    def __init__(self, index=-1, atom1_idx=None, atom2_idx=None,
                 atom3_idx=None, spring_constant=None, eq_angle=None):
        self._index = index
        self._atom1_idx = atom1_idx
        self._atom2_idx = atom2_idx
        self._atom3_idx = atom3_idx
        self._spring_constant = spring_constant
        self._eq_angle = eq_angle

    def set_atom1_idx(self, index):
        self._atom1_idx = index

    def set_atom2_idx(self, index):
        self._atom2_idx = index

    def set_atom3_idx(self, index):
        self._atom3_idx = index

    @property
    def index(self):
        return self._index

    @property
    def atom1_idx(self):
        return self._atom1_idx

    @property
    def atom2_idx(self):
        return self._atom2_idx

    @property
    def atom3_idx(self):
        return self._atom3_idx

    @property
    def spring_constant(self):
        return self._spring_constant

    @property
    def eq_angle(self):
        return self._eq_angle


class Dihedral(TopologyElement):
    _name = 'Dihedral'
    _writable_attrs = ['atom1_idx', 'atom2_idx', 'atom3_idx', 'atom4_idx',
                       'constant', 'prefactor', 'periodicity']

    def __init__(self, index=-1, atom1_idx=None, atom2_idx=None,
                 atom3_idx=None, atom4_idx=None, periodicity=None,
                 prefactor=None, constant=None):
        self._index = index
        self._atom1_idx = atom1_idx
        self._atom2_idx = atom2_idx
        self._atom3_idx = atom3_idx
        self._atom4_idx = atom4_idx
        self._periodicity = periodicity
        self._prefactor = prefactor
        self._constant = constant

    def set_atom1_idx(self, index):
        self._atom1_idx = index

    def set_atom2_idx(self, index):
        self._atom2_idx = index

    def set_atom3_idx(self, index):
        self._atom3_idx = index

    def set_atom4_idx(self, index):
        self._atom4_idx = index

    def plot(self):
        from matplotlib import pyplot
        import numpy as np

        x = unit.Quantity(np.arange(0, np.pi, 0.1), unit=unit.radians)
        pyplot.plot(x, self.constant * (1 + self.prefactor
                                        * np.cos(self.periodicity * x)),
                    'r--')

        pyplot.show()

    @property
    def index(self):
        return self._index

    @property
    def atom1_idx(self):
        return self._atom1_idx

    @property
    def atom2_idx(self):
        return self._atom2_idx

    @property
    def atom3_idx(self):
        return self._atom3_idx

    @property
    def atom4_idx(self):
        return self._atom4_idx

    @property
    def periodicity(self):
        return self._periodicity

    @property
    def prefactor(self):
        return self._prefactor

    @property
    def constant(self):
        return self._constant


class Proper(Dihedral):
    _name = 'Proper'


class Improper(Dihedral):
    _name = 'Improper'


class OFFDihedral(TopologyElement):
    _name = 'OFFDihedral'
    _writable_attrs = ['atom1_idx', 'atom2_idx', 'atom3_idx', 'atom4_idx',
                       'periodicity', 'phase', 'k', 'idivf']
    _to_PELE_class = Dihedral

    def __init__(self, index=-1, atom1_idx=None, atom2_idx=None,
                 atom3_idx=None, atom4_idx=None, periodicity=None,
                 phase=None, k=None, idivf=None):
        self.index = index
        self.atom1_idx = atom1_idx
        self.atom2_idx = atom2_idx
        self.atom3_idx = atom3_idx
        self.atom4_idx = atom4_idx
        self.periodicity = periodicity
        self.phase = phase
        self.k = k
        self.idivf = idivf

    def to_PELE(self):
        if (self.periodicity is None or self.phase is None
                or self.k is None or self.idivf is None):
            return None

        assert self.periodicity in (1, 2, 3, 4, 6), 'Expected values for ' \
            'periodicity are 1, 2, 3, 4 or 6, obtained ' \
            '{}'.format(self.periodicity)
        assert self.phase.value_in_unit(unit.degree) in (0, 180), \
            'Expected values for phase are 0 or 180, obtained ' \
            '{}'.format(self.phase)
        assert self.idivf == 1, 'Expected values for idivf is 1, ' \
            'obtained {}'.format(self.divf)

        if self.phase.value_in_unit(unit.degree) == 180:
            PELE_prefactor = -1
        else:
            PELE_prefactor = 1

        # TODO doublecheck idivf term in OFF's torsion equation
        PELE_constant = self.k / self.idivf

        PELE_dihedral_kwargs = {'index': self.index,
                                'atom1_idx': self.atom1_idx,
                                'atom2_idx': self.atom2_idx,
                                'atom3_idx': self.atom3_idx,
                                'atom4_idx': self.atom4_idx,
                                'periodicity': self.periodicity,
                                'prefactor': PELE_prefactor,
                                'constant': PELE_constant}

        return self._to_PELE_class(**PELE_dihedral_kwargs)

    def plot(self):
        from matplotlib import pyplot
        import numpy as np

        x = unit.Quantity(np.arange(0, np.pi, 0.1), unit=unit.radians)
        pyplot.plot(x,
                    self.k * (1 + np.cos(self.periodicity * x - self.phase)),
                    'r--')

        pyplot.show()


class OFFProper(OFFDihedral):
    _name = 'OFFProper'
    _to_PELE_class = Proper


class OFFImproper(OFFDihedral):
    _name = 'OFFImproper'
    _to_PELE_class = Improper
