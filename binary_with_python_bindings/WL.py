# This file was automatically generated by SWIG (https://www.swig.org).
# Version 4.1.1
#
# Do not make changes to this file unless you know what you are doing - modify
# the SWIG interface file instead.

from sys import version_info as _swig_python_version_info
# Import the low-level C/C++ module
if __package__ or "." in __name__:
    from . import _WL
else:
    import _WL

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)


def _swig_setattr_nondynamic_instance_variable(set):
    def set_instance_attr(self, name, value):
        if name == "this":
            set(self, name, value)
        elif name == "thisown":
            self.this.own(value)
        elif hasattr(self, name) and isinstance(getattr(type(self), name), property):
            set(self, name, value)
        else:
            raise AttributeError("You cannot add instance attributes to %s" % self)
    return set_instance_attr


def _swig_setattr_nondynamic_class_variable(set):
    def set_class_attr(cls, name, value):
        if hasattr(cls, name) and not isinstance(getattr(cls, name), property):
            set(cls, name, value)
        else:
            raise AttributeError("You cannot add class attributes to %s" % cls)
    return set_class_attr


def _swig_add_metaclass(metaclass):
    """Class decorator for adding a metaclass to a SWIG wrapped class - a slimmed down version of six.add_metaclass"""
    def wrapper(cls):
        return metaclass(cls.__name__, cls.__bases__, cls.__dict__.copy())
    return wrapper


class _SwigNonDynamicMeta(type):
    """Meta class to enforce nondynamic attributes (no new attributes) for a class"""
    __setattr__ = _swig_setattr_nondynamic_class_variable(type.__setattr__)


class SwigPyIterator(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _WL.delete_SwigPyIterator

    def value(self):
        return _WL.SwigPyIterator_value(self)

    def incr(self, n=1):
        return _WL.SwigPyIterator_incr(self, n)

    def decr(self, n=1):
        return _WL.SwigPyIterator_decr(self, n)

    def distance(self, x):
        return _WL.SwigPyIterator_distance(self, x)

    def equal(self, x):
        return _WL.SwigPyIterator_equal(self, x)

    def copy(self):
        return _WL.SwigPyIterator_copy(self)

    def next(self):
        return _WL.SwigPyIterator_next(self)

    def __next__(self):
        return _WL.SwigPyIterator___next__(self)

    def previous(self):
        return _WL.SwigPyIterator_previous(self)

    def advance(self, n):
        return _WL.SwigPyIterator_advance(self, n)

    def __eq__(self, x):
        return _WL.SwigPyIterator___eq__(self, x)

    def __ne__(self, x):
        return _WL.SwigPyIterator___ne__(self, x)

    def __iadd__(self, n):
        return _WL.SwigPyIterator___iadd__(self, n)

    def __isub__(self, n):
        return _WL.SwigPyIterator___isub__(self, n)

    def __add__(self, n):
        return _WL.SwigPyIterator___add__(self, n)

    def __sub__(self, *args):
        return _WL.SwigPyIterator___sub__(self, *args)
    def __iter__(self):
        return self

# Register SwigPyIterator in _WL:
_WL.SwigPyIterator_swigregister(SwigPyIterator)
class types_t(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def iterator(self):
        return _WL.types_t_iterator(self)
    def __iter__(self):
        return self.iterator()

    def __nonzero__(self):
        return _WL.types_t___nonzero__(self)

    def __bool__(self):
        return _WL.types_t___bool__(self)

    def __len__(self):
        return _WL.types_t___len__(self)

    def __getslice__(self, i, j):
        return _WL.types_t___getslice__(self, i, j)

    def __setslice__(self, *args):
        return _WL.types_t___setslice__(self, *args)

    def __delslice__(self, i, j):
        return _WL.types_t___delslice__(self, i, j)

    def __delitem__(self, *args):
        return _WL.types_t___delitem__(self, *args)

    def __getitem__(self, *args):
        return _WL.types_t___getitem__(self, *args)

    def __setitem__(self, *args):
        return _WL.types_t___setitem__(self, *args)

    def pop(self):
        return _WL.types_t_pop(self)

    def append(self, x):
        return _WL.types_t_append(self, x)

    def empty(self):
        return _WL.types_t_empty(self)

    def size(self):
        return _WL.types_t_size(self)

    def swap(self, v):
        return _WL.types_t_swap(self, v)

    def begin(self):
        return _WL.types_t_begin(self)

    def end(self):
        return _WL.types_t_end(self)

    def rbegin(self):
        return _WL.types_t_rbegin(self)

    def rend(self):
        return _WL.types_t_rend(self)

    def clear(self):
        return _WL.types_t_clear(self)

    def get_allocator(self):
        return _WL.types_t_get_allocator(self)

    def pop_back(self):
        return _WL.types_t_pop_back(self)

    def erase(self, *args):
        return _WL.types_t_erase(self, *args)

    def __init__(self, *args):
        _WL.types_t_swiginit(self, _WL.new_types_t(*args))

    def push_back(self, x):
        return _WL.types_t_push_back(self, x)

    def front(self):
        return _WL.types_t_front(self)

    def back(self):
        return _WL.types_t_back(self)

    def assign(self, n, x):
        return _WL.types_t_assign(self, n, x)

    def resize(self, *args):
        return _WL.types_t_resize(self, *args)

    def insert(self, *args):
        return _WL.types_t_insert(self, *args)

    def reserve(self, n):
        return _WL.types_t_reserve(self, n)

    def capacity(self):
        return _WL.types_t_capacity(self)
    __swig_destroy__ = _WL.delete_types_t

# Register types_t in _WL:
_WL.types_t_swigregister(types_t)
class coord_t(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def iterator(self):
        return _WL.coord_t_iterator(self)
    def __iter__(self):
        return self.iterator()

    def __nonzero__(self):
        return _WL.coord_t___nonzero__(self)

    def __bool__(self):
        return _WL.coord_t___bool__(self)

    def __len__(self):
        return _WL.coord_t___len__(self)

    def __getslice__(self, i, j):
        return _WL.coord_t___getslice__(self, i, j)

    def __setslice__(self, *args):
        return _WL.coord_t___setslice__(self, *args)

    def __delslice__(self, i, j):
        return _WL.coord_t___delslice__(self, i, j)

    def __delitem__(self, *args):
        return _WL.coord_t___delitem__(self, *args)

    def __getitem__(self, *args):
        return _WL.coord_t___getitem__(self, *args)

    def __setitem__(self, *args):
        return _WL.coord_t___setitem__(self, *args)

    def pop(self):
        return _WL.coord_t_pop(self)

    def append(self, x):
        return _WL.coord_t_append(self, x)

    def empty(self):
        return _WL.coord_t_empty(self)

    def size(self):
        return _WL.coord_t_size(self)

    def swap(self, v):
        return _WL.coord_t_swap(self, v)

    def begin(self):
        return _WL.coord_t_begin(self)

    def end(self):
        return _WL.coord_t_end(self)

    def rbegin(self):
        return _WL.coord_t_rbegin(self)

    def rend(self):
        return _WL.coord_t_rend(self)

    def clear(self):
        return _WL.coord_t_clear(self)

    def get_allocator(self):
        return _WL.coord_t_get_allocator(self)

    def pop_back(self):
        return _WL.coord_t_pop_back(self)

    def erase(self, *args):
        return _WL.coord_t_erase(self, *args)

    def __init__(self, *args):
        _WL.coord_t_swiginit(self, _WL.new_coord_t(*args))

    def push_back(self, x):
        return _WL.coord_t_push_back(self, x)

    def front(self):
        return _WL.coord_t_front(self)

    def back(self):
        return _WL.coord_t_back(self)

    def assign(self, n, x):
        return _WL.coord_t_assign(self, n, x)

    def resize(self, *args):
        return _WL.coord_t_resize(self, *args)

    def insert(self, *args):
        return _WL.coord_t_insert(self, *args)

    def reserve(self, n):
        return _WL.coord_t_reserve(self, n)

    def capacity(self):
        return _WL.coord_t_capacity(self)
    __swig_destroy__ = _WL.delete_coord_t

# Register coord_t in _WL:
_WL.coord_t_swigregister(coord_t)
class coordlist_t(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr

    def iterator(self):
        return _WL.coordlist_t_iterator(self)
    def __iter__(self):
        return self.iterator()

    def __nonzero__(self):
        return _WL.coordlist_t___nonzero__(self)

    def __bool__(self):
        return _WL.coordlist_t___bool__(self)

    def __len__(self):
        return _WL.coordlist_t___len__(self)

    def __getslice__(self, i, j):
        return _WL.coordlist_t___getslice__(self, i, j)

    def __setslice__(self, *args):
        return _WL.coordlist_t___setslice__(self, *args)

    def __delslice__(self, i, j):
        return _WL.coordlist_t___delslice__(self, i, j)

    def __delitem__(self, *args):
        return _WL.coordlist_t___delitem__(self, *args)

    def __getitem__(self, *args):
        return _WL.coordlist_t___getitem__(self, *args)

    def __setitem__(self, *args):
        return _WL.coordlist_t___setitem__(self, *args)

    def pop(self):
        return _WL.coordlist_t_pop(self)

    def append(self, x):
        return _WL.coordlist_t_append(self, x)

    def empty(self):
        return _WL.coordlist_t_empty(self)

    def size(self):
        return _WL.coordlist_t_size(self)

    def swap(self, v):
        return _WL.coordlist_t_swap(self, v)

    def begin(self):
        return _WL.coordlist_t_begin(self)

    def end(self):
        return _WL.coordlist_t_end(self)

    def rbegin(self):
        return _WL.coordlist_t_rbegin(self)

    def rend(self):
        return _WL.coordlist_t_rend(self)

    def clear(self):
        return _WL.coordlist_t_clear(self)

    def get_allocator(self):
        return _WL.coordlist_t_get_allocator(self)

    def pop_back(self):
        return _WL.coordlist_t_pop_back(self)

    def erase(self, *args):
        return _WL.coordlist_t_erase(self, *args)

    def __init__(self, *args):
        _WL.coordlist_t_swiginit(self, _WL.new_coordlist_t(*args))

    def push_back(self, x):
        return _WL.coordlist_t_push_back(self, x)

    def front(self):
        return _WL.coordlist_t_front(self)

    def back(self):
        return _WL.coordlist_t_back(self)

    def assign(self, n, x):
        return _WL.coordlist_t_assign(self, n, x)

    def resize(self, *args):
        return _WL.coordlist_t_resize(self, *args)

    def insert(self, *args):
        return _WL.coordlist_t_insert(self, *args)

    def reserve(self, n):
        return _WL.coordlist_t_reserve(self, n)

    def capacity(self):
        return _WL.coordlist_t_capacity(self)
    __swig_destroy__ = _WL.delete_coordlist_t

# Register coordlist_t in _WL:
_WL.coordlist_t_swigregister(coordlist_t)

def average_stdev(data_list):
    return _WL.average_stdev(data_list)
class configuration(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    N_particles = property(_WL.configuration_N_particles_get, _WL.configuration_N_particles_set)
    v_xyz = property(_WL.configuration_v_xyz_get, _WL.configuration_v_xyz_set)
    v_type = property(_WL.configuration_v_type_get, _WL.configuration_v_type_set)

    def init(self, N_particles_t):
        return _WL.configuration_init(self, N_particles_t)

    def print_config(self):
        return _WL.configuration_print_config(self)

    def perturb(self):
        return _WL.configuration_perturb(self)

    def __init__(self):
        _WL.configuration_swiginit(self, _WL.new_configuration())
    __swig_destroy__ = _WL.delete_configuration

# Register configuration in _WL:
_WL.configuration_swigregister(configuration)

def init_configuration(N_particles_t):
    return _WL.init_configuration(N_particles_t)

def set_configuration(v_xyz_t, v_type_t):
    return _WL.set_configuration(v_xyz_t, v_type_t)

def get_configuration(v_xyz_t, v_type_t):
    return _WL.get_configuration(v_xyz_t, v_type_t)

def rand_gaussian():
    return _WL.rand_gaussian()

def rand_double():
    return _WL.rand_double()
class neighbor(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    member = property(_WL.neighbor_member_get, _WL.neighbor_member_set)
    x_old = property(_WL.neighbor_x_old_get, _WL.neighbor_x_old_set)
    dx = property(_WL.neighbor_dx_get, _WL.neighbor_dx_set)

    def __init__(self):
        _WL.neighbor_swiginit(self, _WL.new_neighbor())
    __swig_destroy__ = _WL.delete_neighbor

# Register neighbor in _WL:
_WL.neighbor_swigregister(neighbor)

def nsq_neighbor_init(x, nbr, cutoff, skin, L):
    return _WL.nsq_neighbor_init(x, nbr, cutoff, skin, L)

def nsq_neighbor_check_fast(x, nbr, cutoff, skin, L):
    return _WL.nsq_neighbor_check_fast(x, nbr, cutoff, skin, L)

def nsq_neighbor_check(x, nbr, cutoff, skin, L):
    return _WL.nsq_neighbor_check(x, nbr, cutoff, skin, L)

def nsq_neighbor_rebuild(x, nbr, cutoff, skin, L):
    return _WL.nsq_neighbor_rebuild(x, nbr, cutoff, skin, L)

def init_system(x, N, density, L, types):
    return _WL.init_system(x, N, density, L, types)

def init_system_binary(x, N, n_particles_0, density, L, types):
    return _WL.init_system_binary(x, N, n_particles_0, density, L, types)

def swap_type(p1_type):
    return _WL.swap_type(p1_type)

def LJ_potential(r2, p1_type, p2_type):
    return _WL.LJ_potential(r2, p1_type, p2_type)

def calc_pe_global(x, types, L, cutoff):
    return _WL.calc_pe_global(x, types, L, cutoff)

def calc_pe_brute(x, types, L, cutoff):
    return _WL.calc_pe_brute(x, types, L, cutoff)

def calc_pe(x, types, L, cutoff):
    return _WL.calc_pe(x, types, L, cutoff)

def check_line(line):
    return _WL.check_line(line)

def load_raw(filename, x):
    return _WL.load_raw(filename, x)

def print_xyz_file(filename, x, types, L):
    return _WL.print_xyz_file(filename, x, types, L)

def print_xyz(dataOut, x, types, L):
    return _WL.print_xyz(dataOut, x, types, L)

def load_xyz(filename, xyz, types):
    return _WL.load_xyz(filename, xyz, types)
class metropolis_NVT(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    dx = property(_WL.metropolis_NVT_dx_get, _WL.metropolis_NVT_dx_set)
    dx_max = property(_WL.metropolis_NVT_dx_max_get, _WL.metropolis_NVT_dx_max_set)
    dx_min = property(_WL.metropolis_NVT_dx_min_get, _WL.metropolis_NVT_dx_min_set)
    dx_target_prob = property(_WL.metropolis_NVT_dx_target_prob_get, _WL.metropolis_NVT_dx_target_prob_set)
    translate_particles = property(_WL.metropolis_NVT_translate_particles_get, _WL.metropolis_NVT_translate_particles_set)
    swap = property(_WL.metropolis_NVT_swap_get, _WL.metropolis_NVT_swap_set)
    swap_target_prob = property(_WL.metropolis_NVT_swap_target_prob_get, _WL.metropolis_NVT_swap_target_prob_set)
    swap_particles = property(_WL.metropolis_NVT_swap_particles_get, _WL.metropolis_NVT_swap_particles_set)
    L = property(_WL.metropolis_NVT_L_get, _WL.metropolis_NVT_L_set)
    sseed = property(_WL.metropolis_NVT_sseed_get, _WL.metropolis_NVT_sseed_set)
    T = property(_WL.metropolis_NVT_T_get, _WL.metropolis_NVT_T_set)
    beta = property(_WL.metropolis_NVT_beta_get, _WL.metropolis_NVT_beta_set)
    time_current = property(_WL.metropolis_NVT_time_current_get, _WL.metropolis_NVT_time_current_set)
    potential = property(_WL.metropolis_NVT_potential_get, _WL.metropolis_NVT_potential_set)
    system_initialized = property(_WL.metropolis_NVT_system_initialized_get, _WL.metropolis_NVT_system_initialized_set)
    N_particles = property(_WL.metropolis_NVT_N_particles_get, _WL.metropolis_NVT_N_particles_set)
    x = property(_WL.metropolis_NVT_x_get, _WL.metropolis_NVT_x_set)
    types = property(_WL.metropolis_NVT_types_get, _WL.metropolis_NVT_types_set)
    configuration_initialized = property(_WL.metropolis_NVT_configuration_initialized_get, _WL.metropolis_NVT_configuration_initialized_set)
    skin = property(_WL.metropolis_NVT_skin_get, _WL.metropolis_NVT_skin_set)
    cutoff = property(_WL.metropolis_NVT_cutoff_get, _WL.metropolis_NVT_cutoff_set)
    nbr = property(_WL.metropolis_NVT_nbr_get, _WL.metropolis_NVT_nbr_set)
    nlist_initialized = property(_WL.metropolis_NVT_nlist_initialized_get, _WL.metropolis_NVT_nlist_initialized_set)
    setup_initialized = property(_WL.metropolis_NVT_setup_initialized_get, _WL.metropolis_NVT_setup_initialized_set)
    xyzFilename = property(_WL.metropolis_NVT_xyzFilename_get, _WL.metropolis_NVT_xyzFilename_set)
    thermoFilename = property(_WL.metropolis_NVT_thermoFilename_get, _WL.metropolis_NVT_thermoFilename_set)
    xyzOut_initialized = property(_WL.metropolis_NVT_xyzOut_initialized_get, _WL.metropolis_NVT_xyzOut_initialized_set)
    thermoOut_initialized = property(_WL.metropolis_NVT_thermoOut_initialized_get, _WL.metropolis_NVT_thermoOut_initialized_set)

    def __init__(self):
        _WL.metropolis_NVT_swiginit(self, _WL.new_metropolis_NVT())

    def init_translate(self, dx_t, dx_min_t, dx_max_t, dx_target_prob_t):
        return _WL.metropolis_NVT_init_translate(self, dx_t, dx_min_t, dx_max_t, dx_target_prob_t)

    def init_system(self, L_t, sseed_t, T_t):
        return _WL.metropolis_NVT_init_system(self, L_t, sseed_t, T_t)

    def init_swap(self, swap_t, swap_target_prob_t):
        return _WL.metropolis_NVT_init_swap(self, swap_t, swap_target_prob_t)

    def init_nlist(self, cutoff_t, skin_t):
        return _WL.metropolis_NVT_init_nlist(self, cutoff_t, skin_t)

    def set_configuration(self, x_t, types_array_t):
        return _WL.metropolis_NVT_set_configuration(self, x_t, types_array_t)

    def get_configuration(self, x_t, types_array_t):
        return _WL.metropolis_NVT_get_configuration(self, x_t, types_array_t)

    def set_traj_filename(self, filename):
        return _WL.metropolis_NVT_set_traj_filename(self, filename)

    def set_thermo_filename(self, filename):
        return _WL.metropolis_NVT_set_thermo_filename(self, filename)

    def set_temperature(self, T_t):
        return _WL.metropolis_NVT_set_temperature(self, T_t)

    def setup(self):
        return _WL.metropolis_NVT_setup(self)

    def converge_window(self, time_run, time_U_update, time_xyz_output, time_thermo_output, adjust, U_min, U_max, T_min, T_max):
        return _WL.metropolis_NVT_converge_window(self, time_run, time_U_update, time_xyz_output, time_thermo_output, adjust, U_min, U_max, T_min, T_max)

    def run(self, time_run, time_U_update, time_xyz_output, time_thermo_output, adjust, append_file):
        return _WL.metropolis_NVT_run(self, time_run, time_U_update, time_xyz_output, time_thermo_output, adjust, append_file)

    def swap_NVT(self, x, types, nbr, L, n_swaps, beta):
        return _WL.metropolis_NVT_swap_NVT(self, x, types, nbr, L, n_swaps, beta)

    def translate_NVT(self):
        return _WL.metropolis_NVT_translate_NVT(self)
    __swig_destroy__ = _WL.delete_metropolis_NVT

# Register metropolis_NVT in _WL:
_WL.metropolis_NVT_swigregister(metropolis_NVT)
class row_container(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    min = property(_WL.row_container_min_get, _WL.row_container_min_set)
    max = property(_WL.row_container_max_get, _WL.row_container_max_set)
    bins = property(_WL.row_container_bins_get, _WL.row_container_bins_set)
    binsize = property(_WL.row_container_binsize_get, _WL.row_container_binsize_set)
    M = property(_WL.row_container_M_get, _WL.row_container_M_set)
    sum = property(_WL.row_container_sum_get, _WL.row_container_sum_set)

    def __init__(self):
        _WL.row_container_swiginit(self, _WL.new_row_container())
    __swig_destroy__ = _WL.delete_row_container

# Register row_container in _WL:
_WL.row_container_swigregister(row_container)
class histogram_2D(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    container = property(_WL.histogram_2D_container_get, _WL.histogram_2D_container_set)
    min_y = property(_WL.histogram_2D_min_y_get, _WL.histogram_2D_min_y_set)
    max_y = property(_WL.histogram_2D_max_y_get, _WL.histogram_2D_max_y_set)
    binsize_y = property(_WL.histogram_2D_binsize_y_get, _WL.histogram_2D_binsize_y_set)
    bins_y = property(_WL.histogram_2D_bins_y_get, _WL.histogram_2D_bins_y_set)
    total_bins = property(_WL.histogram_2D_total_bins_get, _WL.histogram_2D_total_bins_set)
    total_entries = property(_WL.histogram_2D_total_entries_get, _WL.histogram_2D_total_entries_set)
    array = property(_WL.histogram_2D_array_get, _WL.histogram_2D_array_set)
    array_local = property(_WL.histogram_2D_array_local_get, _WL.histogram_2D_array_local_set)
    first_time = property(_WL.histogram_2D_first_time_get, _WL.histogram_2D_first_time_set)
    initialized = property(_WL.histogram_2D_initialized_get, _WL.histogram_2D_initialized_set)

    def __init__(self):
        _WL.histogram_2D_swiginit(self, _WL.new_histogram_2D())

    def init_uniform(self, min_xt, max_xt, binsize_xt, min_yt, max_yt, binsize_yt):
        return _WL.histogram_2D_init_uniform(self, min_xt, max_xt, binsize_xt, min_yt, max_yt, binsize_yt)

    def read_histogram_params(self, filename):
        return _WL.histogram_2D_read_histogram_params(self, filename)

    def check(self):
        return _WL.histogram_2D_check(self)

    def init_y(self, min_yt, max_yt, binsize_yt):
        return _WL.histogram_2D_init_y(self, min_yt, max_yt, binsize_yt)

    def push_back_x(self, min_xt, max_xt, binsize_xt):
        return _WL.histogram_2D_push_back_x(self, min_xt, max_xt, binsize_xt)

    def setup(self):
        return _WL.histogram_2D_setup(self)

    def calculate_min(self):
        return _WL.histogram_2D_calculate_min(self)

    def insert(self, *args):
        return _WL.histogram_2D_insert(self, *args)

    def splat(self, value):
        return _WL.histogram_2D_splat(self, value)

    def print_histogram(self, filename):
        return _WL.histogram_2D_print_histogram(self, filename)

    def print_params(self, filename):
        return _WL.histogram_2D_print_params(self, filename)

    def read_from_file(self, filename, n_threads):
        return _WL.histogram_2D_read_from_file(self, filename, n_threads)

    def clear(self):
        return _WL.histogram_2D_clear(self)

    def clear_all(self):
        return _WL.histogram_2D_clear_all(self)

    def reset(self):
        return _WL.histogram_2D_reset(self)

    def get_hist(self, pp_x, pp_y):
        return _WL.histogram_2D_get_hist(self, pp_x, pp_y)

    def get_index(self, pp_x, pp_y):
        return _WL.histogram_2D_get_index(self, pp_x, pp_y)

    def size(self):
        return _WL.histogram_2D_size(self)

    def n_entries(self):
        return _WL.histogram_2D_n_entries(self)

    def get_min(self, pp_y):
        return _WL.histogram_2D_get_min(self, pp_y)

    def get_max(self, pp_y):
        return _WL.histogram_2D_get_max(self, pp_y)
    __swig_destroy__ = _WL.delete_histogram_2D

# Register histogram_2D in _WL:
_WL.histogram_2D_swigregister(histogram_2D)

def check_flatness(hist, threshold, metric):
    return _WL.check_flatness(hist, threshold, metric)
class WL_2D(object):
    thisown = property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc="The membership flag")
    __repr__ = _swig_repr
    dx = property(_WL.WL_2D_dx_get, _WL.WL_2D_dx_set)
    dx_max = property(_WL.WL_2D_dx_max_get, _WL.WL_2D_dx_max_set)
    dx_min = property(_WL.WL_2D_dx_min_get, _WL.WL_2D_dx_min_set)
    dx_target_prob = property(_WL.WL_2D_dx_target_prob_get, _WL.WL_2D_dx_target_prob_set)
    translate_particles = property(_WL.WL_2D_translate_particles_get, _WL.WL_2D_translate_particles_set)
    nn = property(_WL.WL_2D_nn_get, _WL.WL_2D_nn_set)
    swap = property(_WL.WL_2D_swap_get, _WL.WL_2D_swap_set)
    swap_max = property(_WL.WL_2D_swap_max_get, _WL.WL_2D_swap_max_set)
    swap_target_prob = property(_WL.WL_2D_swap_target_prob_get, _WL.WL_2D_swap_target_prob_set)
    swap_particles = property(_WL.WL_2D_swap_particles_get, _WL.WL_2D_swap_particles_set)
    cswap = property(_WL.WL_2D_cswap_get, _WL.WL_2D_cswap_set)
    cswap_max = property(_WL.WL_2D_cswap_max_get, _WL.WL_2D_cswap_max_set)
    cswap_target_prob = property(_WL.WL_2D_cswap_target_prob_get, _WL.WL_2D_cswap_target_prob_set)
    cswap_particles = property(_WL.WL_2D_cswap_particles_get, _WL.WL_2D_cswap_particles_set)
    vol = property(_WL.WL_2D_vol_get, _WL.WL_2D_vol_set)
    vol_max = property(_WL.WL_2D_vol_max_get, _WL.WL_2D_vol_max_set)
    ln_vol_change_max = property(_WL.WL_2D_ln_vol_change_max_get, _WL.WL_2D_ln_vol_change_max_set)
    vol_target_prob = property(_WL.WL_2D_vol_target_prob_get, _WL.WL_2D_vol_target_prob_set)
    vol_change = property(_WL.WL_2D_vol_change_get, _WL.WL_2D_vol_change_set)
    L = property(_WL.WL_2D_L_get, _WL.WL_2D_L_set)
    sseed = property(_WL.WL_2D_sseed_get, _WL.WL_2D_sseed_set)
    time_current = property(_WL.WL_2D_time_current_get, _WL.WL_2D_time_current_set)
    rank = property(_WL.WL_2D_rank_get, _WL.WL_2D_rank_set)
    threads = property(_WL.WL_2D_threads_get, _WL.WL_2D_threads_set)
    potential = property(_WL.WL_2D_potential_get, _WL.WL_2D_potential_set)
    system_initialized = property(_WL.WL_2D_system_initialized_get, _WL.WL_2D_system_initialized_set)
    N_particles = property(_WL.WL_2D_N_particles_get, _WL.WL_2D_N_particles_set)
    x = property(_WL.WL_2D_x_get, _WL.WL_2D_x_set)
    types = property(_WL.WL_2D_types_get, _WL.WL_2D_types_set)
    x_init = property(_WL.WL_2D_x_init_get, _WL.WL_2D_x_init_set)
    types_init = property(_WL.WL_2D_types_init_get, _WL.WL_2D_types_init_set)
    L_init = property(_WL.WL_2D_L_init_get, _WL.WL_2D_L_init_set)
    configuration_initialized = property(_WL.WL_2D_configuration_initialized_get, _WL.WL_2D_configuration_initialized_set)
    skin = property(_WL.WL_2D_skin_get, _WL.WL_2D_skin_set)
    cutoff = property(_WL.WL_2D_cutoff_get, _WL.WL_2D_cutoff_set)
    nbr = property(_WL.WL_2D_nbr_get, _WL.WL_2D_nbr_set)
    nlist_initialized = property(_WL.WL_2D_nlist_initialized_get, _WL.WL_2D_nlist_initialized_set)
    M_current = property(_WL.WL_2D_M_current_get, _WL.WL_2D_M_current_set)
    M_min = property(_WL.WL_2D_M_min_get, _WL.WL_2D_M_min_set)
    M_max = property(_WL.WL_2D_M_max_get, _WL.WL_2D_M_max_set)
    M_bin = property(_WL.WL_2D_M_bin_get, _WL.WL_2D_M_bin_set)
    U_min = property(_WL.WL_2D_U_min_get, _WL.WL_2D_U_min_set)
    U_max = property(_WL.WL_2D_U_max_get, _WL.WL_2D_U_max_set)
    U_max_global = property(_WL.WL_2D_U_max_global_get, _WL.WL_2D_U_max_global_set)
    U_bin = property(_WL.WL_2D_U_bin_get, _WL.WL_2D_U_bin_set)
    U_bin_global = property(_WL.WL_2D_U_bin_global_get, _WL.WL_2D_U_bin_global_set)
    flatness = property(_WL.WL_2D_flatness_get, _WL.WL_2D_flatness_set)
    f = property(_WL.WL_2D_f_get, _WL.WL_2D_f_set)
    ln_f = property(_WL.WL_2D_ln_f_get, _WL.WL_2D_ln_f_set)
    visited = property(_WL.WL_2D_visited_get, _WL.WL_2D_visited_set)
    ln_dos = property(_WL.WL_2D_ln_dos_get, _WL.WL_2D_ln_dos_set)
    ln_dos_old = property(_WL.WL_2D_ln_dos_old_get, _WL.WL_2D_ln_dos_old_set)
    WL_initialized = property(_WL.WL_2D_WL_initialized_get, _WL.WL_2D_WL_initialized_set)
    xyzFilename = property(_WL.WL_2D_xyzFilename_get, _WL.WL_2D_xyzFilename_set)
    thermoFilename = property(_WL.WL_2D_thermoFilename_get, _WL.WL_2D_thermoFilename_set)
    loggerFilename = property(_WL.WL_2D_loggerFilename_get, _WL.WL_2D_loggerFilename_set)
    xyzOut_initialized = property(_WL.WL_2D_xyzOut_initialized_get, _WL.WL_2D_xyzOut_initialized_set)
    thermoOut_initialized = property(_WL.WL_2D_thermoOut_initialized_get, _WL.WL_2D_thermoOut_initialized_set)
    loggerOut_initialized = property(_WL.WL_2D_loggerOut_initialized_get, _WL.WL_2D_loggerOut_initialized_set)

    def set_thermo_filename(self, filename):
        return _WL.WL_2D_set_thermo_filename(self, filename)

    def set_traj_filename(self, filename):
        return _WL.WL_2D_set_traj_filename(self, filename)
    setup_initialized = property(_WL.WL_2D_setup_initialized_get, _WL.WL_2D_setup_initialized_set)
    n_accept_translate = property(_WL.WL_2D_n_accept_translate_get, _WL.WL_2D_n_accept_translate_set)
    n_accept_swap = property(_WL.WL_2D_n_accept_swap_get, _WL.WL_2D_n_accept_swap_set)
    n_accept_cswap = property(_WL.WL_2D_n_accept_cswap_get, _WL.WL_2D_n_accept_cswap_set)
    n_accept_vol = property(_WL.WL_2D_n_accept_vol_get, _WL.WL_2D_n_accept_vol_set)
    iterations = property(_WL.WL_2D_iterations_get, _WL.WL_2D_iterations_set)
    since_last = property(_WL.WL_2D_since_last_get, _WL.WL_2D_since_last_set)

    def __init__(self):
        _WL.WL_2D_swiginit(self, _WL.new_WL_2D())

    def init_uniform_WL(self, U_min_t, U_max_t, U_bin_t, flatness_t, M_min_t, M_max_t, M_bin_t):
        return _WL.WL_2D_init_uniform_WL(self, U_min_t, U_max_t, U_bin_t, flatness_t, M_min_t, M_max_t, M_bin_t)

    def init_from_params_WL(self, filename, flatness_t):
        return _WL.WL_2D_init_from_params_WL(self, filename, flatness_t)

    def init_ln_dos_from_file(self, filename_ln_dos):
        return _WL.WL_2D_init_ln_dos_from_file(self, filename_ln_dos)

    def init_y(self, M_min_t, M_max_t, M_bin_t):
        return _WL.WL_2D_init_y(self, M_min_t, M_max_t, M_bin_t)

    def push_back_x(self, U_min_t, U_max_t, U_bin_t):
        return _WL.WL_2D_push_back_x(self, U_min_t, U_max_t, U_bin_t)

    def setup_histograms(self, flatness_t):
        return _WL.WL_2D_setup_histograms(self, flatness_t)

    def init_translate(self, dx_t, dx_min_t, dx_max_t, dx_target_prob_t):
        return _WL.WL_2D_init_translate(self, dx_t, dx_min_t, dx_max_t, dx_target_prob_t)

    def init_system(self, L_t, sseed_t):
        return _WL.WL_2D_init_system(self, L_t, sseed_t)

    def init_swap(self, swap_t, swap_max_t, swap_target_prob_t):
        return _WL.WL_2D_init_swap(self, swap_t, swap_max_t, swap_target_prob_t)

    def init_cswap(self, cswap_t, cswap_max_t, cswap_target_prob_t):
        return _WL.WL_2D_init_cswap(self, cswap_t, cswap_max_t, cswap_target_prob_t)

    def init_vol_change(self, vol_t, vol_max_t, vol_target_prob_t, ln_vol_change_max_t):
        return _WL.WL_2D_init_vol_change(self, vol_t, vol_max_t, vol_target_prob_t, ln_vol_change_max_t)

    def init_nlist(self, cutoff_t, skin_t):
        return _WL.WL_2D_init_nlist(self, cutoff_t, skin_t)

    def set_configuration(self, x_t, types_array_t):
        return _WL.WL_2D_set_configuration(self, x_t, types_array_t)

    def reset_configuration(self):
        return _WL.WL_2D_reset_configuration(self)

    def get_configuration(self, x_t, types_array_t):
        return _WL.WL_2D_get_configuration(self, x_t, types_array_t)

    def setup(self):
        return _WL.WL_2D_setup(self)

    def calc_U_min(self, f_t, run_time, time_U_update, time_check_flat, time_thermo_output, adjust, append_file):
        return _WL.WL_2D_calc_U_min(self, f_t, run_time, time_U_update, time_check_flat, time_thermo_output, adjust, append_file)

    def run_step(self, f_t, adjust, append_file, check_now):
        return _WL.WL_2D_run_step(self, f_t, adjust, append_file, check_now)

    def run(self, f_t, time_U_update, time_check_flat, time_xyz_output, time_thermo_output, adjust, append_file):
        return _WL.WL_2D_run(self, f_t, time_U_update, time_check_flat, time_xyz_output, time_thermo_output, adjust, append_file)

    def clear_visited(self):
        return _WL.WL_2D_clear_visited(self)

    def translate_WL(self):
        return _WL.WL_2D_translate_WL(self)

    def swap_WL(self):
        return _WL.WL_2D_swap_WL(self)

    def count_type0(self):
        return _WL.WL_2D_count_type0(self)

    def cswap_WL(self):
        return _WL.WL_2D_cswap_WL(self)

    def calc_vol(self):
        return _WL.WL_2D_calc_vol(self)

    def vol_WL(self):
        return _WL.WL_2D_vol_WL(self)
    __swig_destroy__ = _WL.delete_WL_2D

# Register WL_2D in _WL:
_WL.WL_2D_swigregister(WL_2D)
