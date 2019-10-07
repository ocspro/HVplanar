# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 11:26:33 2019

@author: olechrs
"""
from copy import deepcopy
import transformer_inductance as tf_ind
import coupling_capacitance as tf_cap
import electric_field_distribution as tf_field


def run_FEMM_inductance(geometry, currents, inductances):
    '''
    Launch FEMM to calculate inductance and resistance parameters of a planar
    transformer at 6.78 MHz operating frequency.
    'geometry' - TransformerGeometry object
    'currents' - specification of primary and secondary side currents [i1, i2]
    'inductances' - specification of nono-zero primary and secondary side
    inductances if calculating mutual inductance [L1, L2]
    '''
    geometry.check_variables()
    (L, R) = tf_ind.calc_inductance(geometry, currents, inductances)
    return L, R


def run_FEMM_inductance_prim(geometry):
    ''' Calculate the primary side inductance and resistance '''
    L_prim, R_prim = run_FEMM_inductance(geometry, currents=[1, 0],
                                         inductances=[0, 0])
    return L_prim, R_prim


def run_FEMM_inductance_sec(geometry):
    ''' Calculate the secndary side inductance and resistance '''
    L_sec, R_sec = run_FEMM_inductance(geometry, currents=[0, 1],
                                       inductances=[0, 0])
    return L_sec, R_sec


def run_FEMM_inductance_mutual(geometry):
    ''' Calculate mutual inductance of the transformer at 6.78 MHz operating
    frequency.. First, the primary and secondary self inductances are
    calculated, and the mutual inductance is then found. '''
    L_prim, R_prim = run_FEMM_inductance_prim(geometry)
    L_sec, R_sec = run_FEMM_inductance_sec(geometry)
    L_mutual, __ = run_FEMM_inductance(geometry, currents=[1, 1],
                                       inductances=[L_prim, L_sec])
    return [L_prim, R_prim, L_sec, R_sec, L_mutual]


def run_FEMM_capacitance(geometry):
    '''
    Calculate the coupling capacitance between the primary and secondary
    windings of a planar transformer with geometry of type
    TransformerGeometry
    'geometry' - TransformerGeometry object
    '''
    geometry.check_variables()
    C = tf_cap.calc_coupling_capacitance(geometry)
    return C


def run_FEMM_field_distribution(geometry, voltage=1, view=None):
    '''
    Calculate the field distribution between the primary and secondary
    windings of a planar transformer with geometry of type
    TransformerGeometry
    'geometry' - TransformerGeometry object
    'voltage' - the voltage value of the positive electrode. The negative
    electrode is set to 0. Default value is 1 volt
    'view' - View object
    '''
    geometry.check_variables()
    C = tf_field.calc_field_distribution(geometry, voltage, view)
    return C


class TransformerGeometry(object):
    '''Parameters for defining a planar transformer printed on PCBs with
    an insulating layer sandwiched inbetween.'''
    def __init__(self, turns_primary=1, turns_secondary=1,
                 width_copper=1, width_between_tracks=1, radius_inner_track=4,
                 layers_primary=1, layers_secondary=1, height_pcb_core=.9,
                 height_pcb_prepreg=0, height_dielectric=1,
                 height_copper=0.035, height_gap=0.01, height_gel=5,
                 radius_pcb=None, radius_gel=None, radius_dielectric=None,
                 material_dielectric='fr4', edge_type='Round', guard=None):
        self.turns_primary = turns_primary
        self.layers_primary = layers_primary
        self.turns_secondary = turns_secondary
        self.layers_secondary = layers_secondary
        self.height_pcb_core = height_pcb_core
        self.height_pcb_prepreg = height_pcb_prepreg
        self.height_dielectric = height_dielectric
        self.height_copper = height_copper
        self.height_gap = height_gap
        self.height_gel = height_gel
        self.width_copper = width_copper
        self.width_between_tracks = width_between_tracks
        self.radius_inner_track = radius_inner_track
        self.radius_pcb = radius_pcb
        self.radius_dielectric = radius_dielectric
        self.radius_gel = radius_gel
        self.material_dielectric = material_dielectric
        self.edge_type = edge_type
        self.guard = guard

    def check_variables(self):
        '''
        Check and append changes to the class variables to be consistent
        with the expected format of the methods that use the
        TransformerGeometry class
        '''
        if type(self.width_between_tracks) not in [list, tuple]:
            self.width_between_tracks = [self.width_between_tracks, ] * 2
        if type(self.radius_inner_track) not in [list, tuple]:
            self.radius_inner_track = [self.radius_inner_track, ] * 2
        if self.radius_pcb is None:
            self.radius_pcb = (max(self.turns_primary, self.turns_secondary) *
                               (self.width_between_tracks[0] +
                                self.width_copper) +
                               self.radius_inner_track[0] + 0.5)
        if self.radius_dielectric is None:
            self.radius_dielectric = self.radius_pcb + 20
        if self.radius_gel is None:
            self.radius_gel = self.radius_dielectric + 2


class Guard(object):
    def __init__(self, gtype, polarity):
        self.gtype = gtype
        self.polarity = polarity  # 0 for low, >0 for high

    def __len__(self):
        return 1

    def __str__(self):
        return '<Guard type: {}>'.format(self.gtype)

    def __deepcopy__(self, memo):  # memo is a dict of id's to copies
        id_self = id(self)         # memoization avoids unnecesary recursion
        _copy = memo.get(id_self)
        if _copy is None:
            _copy = type(self)(
                deepcopy(self.gtype, memo),
                deepcopy(self.polarity, memo))
            memo[id_self] = _copy
        return _copy

    def mirror(self):
        _copy = deepcopy(self)
        _copy.polarity = int(not _copy.polarity)
        return _copy

    def pair(self):
        return [self, self.mirror()]


class GuardRing(Guard):
    ''' Parameters for the guard ring structure with a circular cross
    section. '''
    def __init__(self, radius, distance, gap, polarity=1):
        Guard.__init__(self, 'Ring', polarity)
        self.radius = radius
        self.distance = distance
        self.gap = gap

    def __deepcopy__(self, memo):  # memo is a dict of id's to copies
        id_self = id(self)         # memoization avoids unnecesary recursion
        _copy = memo.get(id_self)
        if _copy is None:
            _copy = type(self)(
                deepcopy(self.radius, memo),
                deepcopy(self.distance, memo),
                deepcopy(self.gap, memo),
                deepcopy(self.polarity, memo))
            memo[id_self] = _copy
        return _copy


class GuardTrench(Guard):
    ''' Parameters for the guard trench structure. Common values for
    drill_angle is 90(flat) and 118 (135 in some cases). '''
    def __init__(self, depth, width, distance, polarity=1, outer_angle=10,
                 drill_angle=90, material='Epoxy', inner=False):
        Guard.__init__(self, 'Trench', polarity)
        self.depth = depth
        self.width = width
        self.distance = distance
        self.polarity = polarity
        self.outer_angle = outer_angle
        self.drill_angle = drill_angle
        self.material = material
        self.inner = inner

    def __deepcopy__(self, memo):  # memo is a dict of id's to copies
        id_self = id(self)         # memoization avoids unnecesary recursion
        _copy = memo.get(id_self)
        if _copy is None:
            _copy = type(self)(
                deepcopy(self.depth, memo),
                deepcopy(self.width, memo),
                deepcopy(self.distance, memo),
                deepcopy(self.polarity, memo),
                deepcopy(self.outer_angle, memo),
                deepcopy(self.drill_angle, memo),
                deepcopy(self.material, memo),
                deepcopy(self.inner, memo))
            memo[id_self] = _copy
        return _copy


class View(object):
    ''' Parameters for zoom view in the post-processor of the femm simulation.
    Used for comparable inspection of simulation results across different
    problems and simplfies the work of getting the similar views '''
    filedir = ''  # directory for storing screenshots of View objects

    def __init__(self, span_r, span_z, filename, filetype='bitmap'):
        self.span_r = span_r
        self.span_z = span_z
        self.filename = filename
        self.filetype = filetype
