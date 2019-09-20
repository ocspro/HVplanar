# -*- coding: utf-8 -*-
"""
Created on Tue Apr 03 15:30:59 2018

@author: olechrs
"""


class TransformerGeometry(object):
    ''' Parameters for defining a planar transformer with insulating layer'''

    def __init__(self, turns_primary=1, turns_secondary=1,
                 width_copper=1, width_between_tracks=1, radius_inner_track=4,
                 layers_primary=1, layers_secondary=1, height_pcb_core=.9,
                 height_pcb_prepreg=0, height_dielectric=1,
                 height_copper=0.035, height_gap=0.01, height_gel=5,
                 radius_pcb=None, radius_gel=None, radius_dielectric=100,
                 material_dielectric='fr4'):
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