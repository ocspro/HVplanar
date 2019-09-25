# -*- coding: utf-8 -*-
"""
Created on Tue Apr 03 15:30:59 2018

@author: olechrs
"""

from random import uniform
import transformer_inductance
import coupling_capacitance
import highvoltage_field_distribution as hv_fd
import inspyred


class OptimizationProblem(object):

    def __init__(self):
        """ bla bla """
        self.popsize = 0
        self.maxgeneneration = 0
        self.goal = []
        self.lower_bound = []
        self.upper_bound = []

    def __del__(self):
        pass

    def generator(self, random, args):
        '''
        Initialize and create a random initial population with a uniform
        distribution.
        '''
        chromosome = []
        for lo, hi in zip(self.lower_bound, self.upper_bound):
            chromosome.append(uniform(lo, hi))
        return chromosome

    def evaluator_inductance(self, candidates, args):
        ''' bla bla '''
        fitness = []
        for c in candidates:
            c[:2] = [int(x) for x in c[:2]]
            fit = self.run_FEMM_inductance_prim(c, **args)
            values = []
            for f, g in zip(fit, self.goal):
                if g is None:
                    values.append(f)
                else:
                    values.append(abs(g-f))
            fitness.append(inspyred.ec.emo.Pareto(values))
        return fitness

    def run_FEMM_inductance(self, chromosones, **kwargs):
        '''
        Launch FEMM to calculate inductance and resistance parameters of a
        planar transformer. **kwargs must include
        'geometry' - TransformerGeometry object
        'optimisation_vaiables' - string names of variables in
            TransformerGeometry variable that are given by chromosones in
            order of apperance in chromosones
        'currents' - specification of primary and secondary side currents
            [i1, i2]
        'inductances' - specification of nono-zero primary and secondary side
            inductances if calculating mutual inductance [l1, l2]
        '''
        geometry = kwargs['geometry']
        for var, c in zip(kwargs['optimisation_variables'], chromosones):
            setattr(geometry, var, c)
        if type(geometry.width_between_tracks) not in [list, tuple]:
            geometry.width_between_tracks = \
                [geometry.width_between_tracks, ] * 2
        if type(geometry.radius_inner_track) not in [list, tuple]:
            geometry.radius_inner_track = \
                [geometry.radius_inner_track, ] * 2
        geometry.radius_pcb = (max(geometry.turns_primary,
                                   geometry.turns_secondary) *
                               (geometry.width_between_tracks[0] +
                                geometry.width_copper) +
                               geometry.radius_inner_track[0] + 2)
        if geometry.radius_gel is None:
            geometry.radius_gel = geometry.radius_dielectric + 2
        (L, R) = transformer_inductance.calc_inductance(
            geometry, kwargs['currents'], kwargs['inductances'])
        return L, R

    def run_FEMM_inductance_prim(self, chromosones, **kwargs):
        ''' Calculate the primary side inductance and resistance '''
        l_prim, r_prim = self.run_FEMM_inductance(chromosones, currents=[1, 0],
                                                  inductances=[0, 0], **kwargs)
        return l_prim, r_prim

    def run_FEMM_inductance_sec(self, chromosones, **kwargs):
        ''' Calculate the secndary side inductance and resistance '''
        l_sec, r_sec = self.run_FEMM_inductance(chromosones, currents=[0, 1],
                                                inductances=[0, 0], **kwargs)
        return l_sec, r_sec

    def run_FEMM_inductance_mutual(self, chromosones, **kwargs):
        ''' Calculate mutual inductance of the transformer. First, the primary
        and secondary self inductances are calculated, and the mutual
        inductance is then found.
        '''
        l_prim, r_prim = self.run_FEMM_inductance_prim(chromosones, **kwargs)
        l_sec, r_sec = self.run_FEMM_inductance_sec(chromosones, **kwargs)
        l_mutual, r = self.run_FEMM_inductance(chromosones, currents=[1, 1],
                                               inductances=[l_prim, l_sec],
                                               **kwargs)
        return [l_prim, r_prim, l_sec, r_sec, l_mutual]

    def run_FEMM_capacitance(self, chromosones, **kwargs):
        '''
        Calculate the coupling capacitance between the primary and secondary
        windings of a planar transformer with geometry of type
        TransformerGeometry
        'geometry' - TransformerGeometry object
        'optimisation_vaiables' - string names of variables in
            TransformerGeometry variable that are given by chromosones in
            order of apperance in chromosones
        '''
        geometry = kwargs['geometry']
        for var, c in zip(kwargs['optimisation_variables'], chromosones):
            setattr(geometry, var, c)
        if geometry.radius_pcb is None:
            geometry.radius_pcb = (max(geometry.turns_primary,
                                       geometry.turns_secondary) *
                                   (geometry.width_between_tracks +
                                   geometry.width_copper) +
                                   geometry.radius_inner_track + 2)
        if type(geometry.width_between_tracks) not in [list, tuple]:
            geometry.width_between_tracks = \
                [geometry.width_between_tracks, ] * 2
        if type(geometry.radius_inner_track) not in [list, tuple]:
            geometry.radius_inner_track = \
                [geometry.radius_inner_track, ] * 2
        if geometry.radius_gel is None:
            geometry.radius_gel = geometry.radius_dielectric + 2
        cap = coupling_capacitance.calc_coupling_capacitance(geometry)
        return cap

    def run_FEMM_field_distribution(self, chromosones, **kwargs):
        '''
        Calculate the field distribution between the primary and secondary
        windings of a planar transformer with geometry of type
        TransformerGeometry
        'geometry' - TransformerGeometry object
        'optimisation_vaiables' - string names of variables in
            TransformerGeometry variable that are given by chromosones in
            order of apperance in chromosones
        'voltage' - the voltage value of the positive electrode. The negative
        electrode is set to 0. Default value is 1 volt
        '''
        voltage = 1
        filedir = None
        guard = None
        geometry = kwargs['geometry']
        for var, c in zip(kwargs['optimisation_variables'], chromosones):
            setattr(geometry, var, c)
        if geometry.radius_pcb is None:
            geometry.radius_pcb = (max(geometry.turns_primary,
                                       geometry.turns_secondary) *
                                   (geometry.width_between_tracks +
                                   geometry.width_copper) +
                                   geometry.radius_inner_track + .25)
        if type(geometry.width_between_tracks) not in [list, tuple]:
            geometry.width_between_tracks = \
                [geometry.width_between_tracks, ] * 2
        if type(geometry.radius_inner_track) not in [list, tuple]:
            geometry.radius_inner_track = \
                [geometry.radius_inner_track, ] * 2
        if geometry.radius_gel is None:
            geometry.radius_gel = geometry.radius_dielectric + 2
        if 'voltage' in kwargs:
            voltage = kwargs['voltage']
        if 'filedir' in kwargs:
            filedir = kwargs['filedir']
        if 'guard' in kwargs:
            guard = kwargs['guard']
        cap = hv_fd.calc_field_distribution(geometry, voltage, guard, filedir)
        files = [filedir + '\\output_field_edge' + str(idx) + '.csv' for idx
                 in range(4)]
        f_max = [hv_fd.get_field_contour_max(f) for f in files]
        if guard:
            f_max.append(
            hv_fd.get_field_contour_max(filedir + '\\output_field_guard0.csv'))
        else:
            f_max.append(0)

        return cap, f_max
