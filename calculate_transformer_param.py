# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 10:24:20 2018

@author: olechrs
"""

import sys
import numpy as np
import OptimizationProblem as op
import TransformerGeometry as tg


def transformer_inductance(layers, turns, radius_inner, h_dielec, pcb_core,
                           material):
    p = op.OptimizationProblem()
    g = tg.TransformerGeometry(layers_primary=layers,
                               layers_secondary=layers,
                               turns_primary=turns,
                               turns_secondary=turns,
                               height_pcb_core=pcb_core,
                               height_pcb_prepreg=.27,
                               height_dielectric=h_dielec,
                               width_copper=3,
                               width_between_tracks=.5,
                               radius_inner_track=radius_inner,
                               radius_dielectric=radius_inner+50,
                               material_dielectric=material)
    l, r = p.run_FEMM_inductance([], optimisation_variables=[], geometry=g,
                                 currents=[1, 0], inductances=[0, 0])
    print ('L: {:.3f}, R: {:.3f}'.format(l*1e6, r))


def transformer_capacitance(layers, turns, radius_inner, h_dielec,
                            pcb_core, material):
    ''' '''
    p = op.OptimizationProblem()
    g = tg.TransformerGeometry(layers_primary=layers,
                               layers_secondary=layers,
                               turns_primary=turns,
                               turns_secondary=turns,
                               height_pcb_core=pcb_core,
                               height_pcb_prepreg=.27,
                               height_dielectric=h_dielec,
                               width_copper=3,
                               width_between_tracks=.5,
                               radius_inner_track=radius_inner,
                               radius_dielectric=radius_inner+50,
                               material_dielectric=material)
    c = p.run_FEMM_capacitance([], optimisation_variables=[], geometry=g)
    print ('C: {:.1f}'.format(c*1e12))


def transformer_coupling(layers, max_turns, radius_inner, h_dielec, pcb_core,
                         material):
    ''' '''
    p = op.OptimizationProblem()
    g = tg.TransformerGeometry(layers_primary=layers,
                               layers_secondary=layers,
                               turns_primary=turns,
                               turns_secondary=turns,
                               height_pcb_core=pcb_core,
                               height_pcb_prepreg=.27,
                               height_dielectric=h_dielec,
                               width_copper=3,
                               width_between_tracks=.5,
                               radius_inner_track=radius_inner,
                               radius_dielectric=radius_inner+50,
                               material_dielectric=material)
    l1, r1, l2, r2, m = p.run_FEMM_inductance_mutual([],
                                                     optimisation_variables=[],
                                                     geometry=g)
    k = m / np.sqrt(l1 * l2)
    print ('k: {:.3f}'.format(k))


def transformer_all_parameters(layers, turns, radius_inner, h_dielec, pcb_core,
                               material):
    ''' Find all parameters of a single transformer design '''
    p = op.OptimizationProblem()

    g = tg.TransformerGeometry(layers_primary=layers,
                               layers_secondary=layers,
                               turns_primary=turns,
                               turns_secondary=turns,
                               height_pcb_core=pcb_core,
                               height_pcb_prepreg=.27,
                               height_dielectric=h_dielec,
                               width_copper=3,
                               width_between_tracks=.5,
                               radius_inner_track=radius_inner,
                               radius_dielectric=radius_inner+50,
                               material_dielectric=material)
    l1, r1, l2, r2, m = p.run_FEMM_inductance_mutual([],
                                                     optimisation_variables=[],
                                                     geometry=g)
    k = m / np.sqrt(l1 * l2)
    c = p.run_FEMM_capacitance([], optimisation_variables=[], geometry=g)
    print ('Lp: {:.3f}, Rp: {:.3f}'.format(l1*1e6, r1))
    print ('Ls: {:.3f}, Rs: {:.3f}'.format(l2*1e6, r2))
    print ('M: {:.3f}, k: {:.3f}'.format(m*1e6, k))
    print ('C: {:.1f}'.format(c*1e12))


def transformer_field(layers, turns, radius_inner, h_dielec,
                      pcb_core, material, gap=None, voltage=60e3, guard=None):
    '''
    Guard ring is defined as follows: (x coordinate, y coordinate, radius)
    Guard ring can be given as a list?
    '''
    p = op.OptimizationProblem()

    g = tg.TransformerGeometry(layers_primary=layers,
                               layers_secondary=layers,
                               turns_primary=turns,
                               turns_secondary=turns,
                               height_pcb_core=pcb_core,
                               height_pcb_prepreg=.24,
                               height_dielectric=h_dielec,
                               width_copper=3,
                               width_between_tracks=.5,
                               radius_inner_track=radius_inner,
                               radius_dielectric=radius_inner+50,
                               material_dielectric=material)
    if gap:
        g.height_gap = gap
    c, f_max = p.run_FEMM_field_distribution([],
                                             optimisation_variables=[],
                                             geometry=g,
                                             voltage=voltage,
                                             filedir=('M:\\PhD\\Paper EPE 19'
                                                      '\\HVplanar'),
                                             guard=guard)
    f_max_list = ', '.join(['{:.1f}'.format(f*1e-6) for f in f_max])
    print ('C: {:.1f}'.format(c*1e12/60e3))
    print ('Peak E-field - upper inner/outer, lower inner/outer, guard:')
    print (f_max_list)

# command structure:
# calculate_transformer_param [layers] [height] [command] [material]
#                             [max_turns] [gap height] [guard]
layers = int(sys.argv[1])
d_height = float(sys.argv[2])
command = int(sys.argv[3])
gap = None
guard = None

if len(sys.argv) > 4:
    material = sys.argv[4].lower()
else:
    material = 'teflon'

turns = int(sys.argv[5])

if len(sys.argv) > 6:
    gap = max(.01, float(sys.argv[6]))

if len(sys.argv) > 7:
    guard = int(sys.argv[7])
    if len(sys.argv) > 8:
        radius = float(sys.argv[8])
    else:
        radius = 1
    if guard == 0:
        guard = [(0.05, 0, radius, 1), (0.05, 0, radius, 0)]
    elif guard == 1:
        guard = [(0.05, -gap, radius, 1), (0.05, gap, radius, 0)]
    elif guard == 2:
        guard = [(0.05, -gap+radius/2., radius, 1), (0.05, gap-radius/2.,
                                                     radius, 0)]

if layers == 1:
    pcb_core_height = .9
    radius_inner = 11.4
    turns = int(sys.argv[5])
elif layers == 2:
    pcb_core_height = .9
    # radii values for other h_dielectric
    radius_inner = 9.2
#    radius_inner = 17.8
#    radius_inner = 49.5

elif layers == 4:
    max_turns = 4
    pcb_core_height = .33
    # radii values for other h_dielectric
#    radius_inner = 16.2
    radius_inner = 5.7

if command == 1:
    transformer_inductance(layers, turns, radius_inner, d_height,
                           pcb_core_height, material)
elif command == 2:
    transformer_capacitance(layers, turns, radius_inner, d_height,
                            pcb_core_height, material)
elif command == 3:
    transformer_coupling(layers, turns, radius_inner, d_height,
                         pcb_core_height, material)
elif command == 4:
    transformer_all_parameters(layers, turns, radius_inner, d_height,
                               pcb_core_height, material)
elif command == 5:
    res = transformer_field(layers, turns, radius_inner, d_height,
                            pcb_core_height, material, gap=gap, guard=guard)
# hv.plot_field_contour(['output_field_edge' + str(idx) +'.csv'
#                        for idx in range(4)])
# hv.plot_field_contour(('output_field_guard0.csv',))
