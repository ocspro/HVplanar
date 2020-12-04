#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module for calculating inducances, resistances and coupling capacitance of
a planar transformer consisting of PCB printed windigns on both sides. Uses the
program FEMM (www.femm.info) to perform finite element calculations. This
module is just simplifying the building of the geometry in FEMM.
"""

from copy import deepcopy

import transformer_inductance as tf_ind
import electric_field_distribution as tf_field


def _run_FEMM_inductance(geometry, currents, inductances, **kwargs):
    '''Launch FEMM to calculate inductance and resistance parameters of a
    planar transformer at 6.78 MHz operating frequency.

    Args:
        geometry (:obj:'TransformerGeometry'): geometry of the planar
        transformer.
        currents (:obj:'list' of int): specification of primary and secondary
        side currents in the form [I_prim, I_sec].
        inductances (:obj:'list' of int): specification of non-zero primary
        and secondary side inductances if calculating mutual inductance
        [L_prim, L_sec].

    Returns:
        inductance (int): calculated inductance value
        resistance (int): calculated resistance value
    '''

    geometry.check_variables()
    (inductance, resistance) = tf_ind.calc_inductance(geometry, currents,
                                                      inductances, **kwargs)
    return inductance, resistance


def run_FEMM_inductance_prim(geometry, **kwargs):
    '''Calculate the primary side inductance and resistance.

    Args:
        geometry (:obj:'TransformerGeometry'): geometry of the planar
        transformer.

    Returns:
        L_prim (int): calculated self-inductance value of primary side
        R_prim (int): calculated eqv. series resistance of primary side
    '''

    L_prim, R_prim = _run_FEMM_inductance(geometry, currents=[1, 0],
                                          inductances=[0, 0], **kwargs)
    return L_prim, R_prim


def run_FEMM_inductance_sec(geometry, **kwargs):
    '''Calculate the secndary side inductance and resistance.

    Args:
        geometry (:obj:'TransformerGeometry'): geometry of the planar
        transformer.

    Returns:
        L_sec (int): calculated self-inductance value of secondary side
        R_sec (int): calculated eqv. series resistance of secondary side
    '''

    L_sec, R_sec = _run_FEMM_inductance(geometry, currents=[0, 1],
                                        inductances=[0, 0], **kwargs)
    return L_sec, R_sec


def run_FEMM_inductance_mutual(geometry, **kwargs):
    '''.Calculate mutual inductance of the transformer.

    Args:
        geometry (:obj:'TransformerGeometry'): geometry of the planar
        transformer.

    Returns:
        L_prim (int): calculated self-inductance value of primary side
        R_prim (int): calculated eqv. series resistance of primary side
        L_sec (int): calculated self-inductance value of secondary side
        R_sec (int): calculated eqv. series resistance of secondary side
        L_mutual (int): calculated  mutual inductance of transformer
    '''

    L_prim, R_prim = run_FEMM_inductance_prim(geometry, **kwargs)
    L_sec, R_sec = run_FEMM_inductance_sec(geometry, **kwargs)
    L_mutual, __ = _run_FEMM_inductance(geometry, currents=[1, 1],
                                        inductances=[L_prim, L_sec], **kwargs)
    return [L_prim, R_prim, L_sec, R_sec, L_mutual]


def run_FEMM_capacitance(geometry, **kwargs):
    '''Calculate the coupling capacitance between the primary and secondary
    windings of the planar transformer.

    Args:
        geometry (:obj:'TransformerGeometry'): geometry of the planar
        transformer.

    Returns:
        capacitance (int): calculated coupling capacitance of transformer.
    '''

    geometry.check_variables()
    capacitance = tf_field.calc_field_distribution(geometry, **kwargs)
    return capacitance


def run_FEMM_field_distribution(geometry, voltage=1, view=None, **kwargs):
    '''
    Calculate the coupling capacitance and investigate the electric field
    distribution between the primary and secondary windings of a planar
    transformer.

    Args:
        geometry (:obj:'TransformerGeometry'): geometry of the planar
        transformer.
        voltage (int): the voltage value of the positive electrode. The
        negative electrode is set to 0.
        view (:obj:'list' of :obj:'View'): View object

    Returns:
        capacitance (int): calculated coupling capacitance of transformer.
    '''

    geometry.check_variables()
    capacitance = tf_field.calc_field_distribution(geometry,
                                                   voltage,
                                                   view,
                                                   **kwargs)
    return capacitance


class TransformerGeometry(object):
    '''Parameters for defining a planar transformer printed on PCBs with
    an insulating layer sandwiched inbetween.'''
    def __init__(self, turns_primary=1, turns_secondary=1,
                 width_copper=1, width_between_tracks=1, radius_inner_track=4,
                 layers_primary=1, layers_secondary=1, height_pcb_core=.9,
                 height_pcb_prepreg=0, height_dielectric=1,
                 height_copper=0.035, height_gap=0.01, height_gel=5,
                 radius_pcb=None, radius_gel=None, radius_dielectric=None,
                 material_dielectric='fr4', tracks=None, guard=None):
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
        # self.edge_type = edge_type
        self.tracks = tracks
        self.guard = guard

    def check_variables(self):
        '''Check and append changes to the class variables to be consistent
        with the expected format of the methods that use this class.'''
        if not isinstance(self.width_between_tracks, (list, tuple)):
            self.width_between_tracks = [self.width_between_tracks, ] * 2
        if not isinstance(self.width_copper, (list, tuple)):
            self.width_between_tracks = [self.width_copper, ] * 2
        if not isinstance(self.radius_inner_track, (list, tuple)):
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


class FancyTrack(object):
    '''Additional manipulation of individual tracks in the winding
    structure.'''
    def __init__(self, layer, track, side_h, side_v='high', rounding='single',
                 elongation=None):
        self.layer = layer
        self.track = track
        self.side_h = side_h
        self.side_v = side_v
        self.rounding = rounding
        self.elongation = elongation

    def __str__(self):
        return ('FancyTrack on {} side, layer {}, track {} with {} rounding.'
                .format(self.side_v, self.layer, self.track, self.rounding))

    def __deepcopy__(self, memo):  # memo is a dict of id's to copies
        id_self = id(self)         # memoization avoids unnecesary recursion
        _copy = memo.get(id_self)
        if _copy is None:
            _copy = type(self)(
                deepcopy(self.layer, memo),
                deepcopy(self.track, memo),
                deepcopy(self.side_h, memo),
                deepcopy(self.side_v, memo),
                deepcopy(self.rounding, memo),
                deepcopy(self.elongation, memo))
            memo[id_self] = _copy
        return _copy

    def mirror(self):
        'Returns a new object with the opposite side than self.'''
        _copy = deepcopy(self)
        if _copy.side_h == 'high':
            _copy.side_h == 'low'
        else:
            _copy.side_h = 'high'
        return _copy


class Guard(object):
    ''' Base class for field guards with different cross section geometry.'''
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
        '''Returns a new guard object with the opposite polarity than self.'''
        _copy = deepcopy(self)
        _copy.polarity = int(not _copy.polarity)
        return _copy

    def pair(self):
        '''Returns a list containing self and a mirrored guard.'''
        return [self, self.mirror()]


class GuardRing(Guard):
    ''' Guard structure with a circular cross section.'''
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
    '''Guard structure with a squarish cross section. Common values for
    drill_angle is 90(flat) and 118 (135 in some cases).'''
    def __init__(self, depth, width, distance, polarity=1, outer_angle=10,
                 drill_angle=90, material='epoxy', inner=False):
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
