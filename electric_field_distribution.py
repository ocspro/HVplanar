#!/usr/bin/env python
# -*- coding: utf-8 -*-

from math import pi
import numpy as np
import matplotlib.pyplot as plt

import femm


def initial_setup(limite_externe, voltage_high, **kwargs):
    '''Start femm and setup problem definition and set boundary condition.'''
    if kwargs.get('hide') is True:
        femm.openfemm(1)
        femm.main_minimize()
    else:
        femm.openfemm()

    # newdocument(doctype)
    # From manual: Creates a new preprocessor document and opens up a new
    # preprocessor window. Specify doctype to be 0 for a magnetics problem, 1
    # for an electrostatics problem, 2 for a heat flow problem, or 3 for a
    # current flow problem. An alternative syntax for this command is
    # create(doctype)
    femm.newdocument(1)

    # ei probdef(units,type,precision,(depth),(minangle))
    # From manual: changes the problem definition. The units parameter
    # specifies the units used for measuring length in the problem domain.
    # Valid "units" entries are "inches", "millimeters", "centimeters", "mils",
    # "meters, and "micrometers". Set problemtype to "planar" for a 2-D planar
    # problem, or to "axi" for an axisymmetric problem. The precision parameter
    # dictates the precision required by the solver. For example, entering
    # 1.E-8 requires the RMS of the residual to be less than 10^-8. A fourth
    # parameter, representing the depth of the problem in the into-thepage
    # direction for 2-D planar problems, can also be specified for planar
    # problems. A sixth parameter represents the minimum angle constraint sent
    # to the mesh generator.
    if 'precision' in kwargs:
        precision = kwargs['precision']
    else:
        precision = 1e-9
    if 'min_angle' in kwargs:
        min_angle = kwargs['min_angle']
    else:
        min_angle = 30
    femm.ei_probdef('millimeters', 'axi', precision, 100, min_angle)

    # Circuit parameters
    # From manual: ei_addconductorprop("conductorname", Vc, qc, conductortype)
    # adds a new conductor property with name "conductorname" with either a
    # prescribed voltage or a prescribed total charge. Set the unused property
    # to zero. The conductortype parameter is 0 for prescribed charge and 1 for
    # prescribed voltage
    femm.ei_addconductorprop('high', voltage_high, 0, 1)
    femm.ei_addconductorprop('zero', 0, 0, 1)

    # Add materials properties used in the simulation
    # ei_addmaterial(’matname’, ex, ey, qv)
    # From manual: adds a new material with called ’matname’ with the material
    # properties:
    # ex Relative permittivity in the x- or r-direction.
    # ey Relative permittivity in the y- or z-direction.
    # qv Volume charge density in units of C/m3
    femm.ei_addmaterial('air', 1, 1, 0)
    femm.ei_addmaterial('fr4', 4.4, 4.4, 0)
    femm.ei_addmaterial('polysterimide', 3.5, 3.5, 0)
    femm.ei_addmaterial('teflon', 2.1, 2.1, 0)
    femm.ei_addmaterial('silgel', 2.7, 2.7, 0)
    femm.ei_addmaterial('midel', 3.15, 3.15, 0)
    if 'material' in kwargs:
        for m in kwargs['material']:
            femm.ei_addmaterial(*m)

    # Boundary conditions
    # ei makeABC(n,R,x,y,bc)
    # From manual: creates a series of circular shells that emulate the
    # impedance of an unbounded domain (i.e. an Improvised Asymptotic Boundary
    # Condition). The n parameter contains the number of shells to be used
    # (should be between 1 and 10), R is the radius of the solution domain, and
    # (x,y) denotes the center of the solution domain. The bc parameter should
    # be specified as 0 for a Dirichlet outer edge or 1 for a Neumann outer
    # edge. If the function is called without all the parameters, the function
    # makes up reasonable values for the missing parameters.
    femm.ei_makeABC(7, limite_externe, 0, 0, 0)


def coords_rectangle(x0, y0, dx, dy):
    '''Return array with the coordinates of the two opposite corners of a
    rectangle.'''
    return (x0, y0, x0+dx, y0+dy)


def get_z_coord_copper_layer(g):
    '''Calculate z coordinate of copper layers depending on the number of
    copper layers in each pcb.'''
    z_coord_first_layer = (g.height_dielectric + g.height_gap, -g.height_gap)
    z_coords = []
    if g.layers_primary == 2:
        z_coords.append((z_coord_first_layer[0],
                         z_coord_first_layer[0] + g.height_copper +
                         g.height_pcb_core))
    elif g.layers_primary == 4:
        z_coords.append((z_coord_first_layer[0],
                         z_coord_first_layer[0] + g.height_copper +
                         g.height_pcb_prepreg,
                         z_coord_first_layer[0] + g.height_copper * 2 +
                         g.height_pcb_prepreg + g.height_pcb_core,
                         z_coord_first_layer[0] + g.height_copper * 3 +
                         g.height_pcb_prepreg * 2 + g.height_pcb_core))
    else:
        z_coords.append((z_coord_first_layer[0],))

    if g.layers_secondary == 2:
        z_coords.append((z_coord_first_layer[1],
                         z_coord_first_layer[1] - g.height_copper -
                         g.height_pcb_core))
    elif g.layers_secondary == 4:
        z_coords.append((z_coord_first_layer[1],
                         z_coord_first_layer[1] - g.height_copper -
                         g.height_pcb_prepreg,
                         z_coord_first_layer[1] - g.height_copper * 2 -
                         g.height_pcb_prepreg - g.height_pcb_core,
                         z_coord_first_layer[1] - g.height_copper * 3 -
                         g.height_pcb_prepreg * 2 - g.height_pcb_core))
    else:
        z_coords.append((z_coord_first_layer[1],))

    return z_coords


def set_conductor_boundry(coords, in_conductor):
    '''Set voltage potential of the four sides of a square conductor.'''
    femm.ei_selectsegment(coords[0], np.average((coords[1], coords[3])))
    femm.ei_selectsegment(np.average((coords[0], coords[2])), coords[1])
    femm.ei_selectsegment(coords[2], np.average((coords[1], coords[3])))
    femm.ei_selectsegment(np.average((coords[0], coords[2])), coords[3])
    if in_conductor == 1:
        femm.ei_setsegmentprop('<None>', 0, 1, 0, 0, 'high')
    elif in_conductor == 0:
        femm.ei_setsegmentprop('<None>', 0, 1, 0, 0, 'zero')
    femm.ei_clearselected()


def set_conductor_label(coords, label_name, label_dict):
    '''Set a block label inside a square conductor and set label properties.'''
    label_coord = (np.average((coords[0], coords[2])),
                   np.average((coords[1], coords[3])))
    femm.ei_addblocklabel(*label_coord)
    femm.ei_selectlabel(*label_coord)
    # No mesh is needed inside conductors
    femm.ei_setblockprop('<No Mesh>', 1, 0, 0)
    femm.ei_clearselected()
    label_dict[label_name] = label_coord


def draw_conductor(coords, in_conductor, label_name, label_dict):
    '''Draw rectangular conductor, set boundary conductions and add label.'''
    femm.ei_drawrectangle(*coords)
    set_conductor_boundry(coords, in_conductor)
    set_conductor_label(coords, label_name, label_dict)


def add_conductors(tg, label_dict):
    '''Add all conductors of primary and secondary side pcbs

    Args:
        tg (TransformerGeometry): contains all geometry information needed
        to build conductors in the problem.
        label_dict (dict): stores label coordinates
    '''

    z_copper = get_z_coord_copper_layer(tg)
    # adding conductors on the primary side
    for c in range(tg.layers_primary):
        for i in range(tg.turns_primary):
            coords_conductor = coords_rectangle(tg.radius_inner_track[0] +
                                                (tg.width_between_tracks[0] +
                                                 tg.width_copper[0]) * i,
                                                z_copper[0][c],
                                                tg.width_copper[0],
                                                tg.height_copper)
            draw_conductor(coords_conductor,
                           1,
                           'prim' + str(i + c * tg.turns_primary),
                           label_dict)
    # adding conductors on the secondary side
    for c in range(tg.layers_secondary):
        for i in range(tg.turns_secondary):
            coords_conductor = coords_rectangle(tg.radius_inner_track[1] +
                                                (tg.width_between_tracks[1] +
                                                 tg.width_copper[1]) * i,
                                                z_copper[1][c],
                                                tg.width_copper[1],
                                                -tg.height_copper)
            draw_conductor(coords_conductor,
                           0,
                           'sec' + str(i + c * tg.turns_secondary),
                           label_dict)


def add_pcbs(tg, label_dict):
    '''Add primary and secondary side pcbs.

    Args:
        tg (TransformerGeometry): contains all geometry information needed
        to build conductors in the problem.
        label_dict (dict): stores coordinates
    '''

    # coordinates of the PCB on the primary side
    if tg.layers_primary == 1 or tg.layers_primary == 2:
        coords_pcb_prim = coords_rectangle(0,
                                           tg.height_dielectric +
                                           tg.height_copper + tg.height_gap,
                                           tg.radius_pcb,
                                           tg.height_pcb_core)
    elif tg.layers_primary == 4:
        coords_pcb_prim = coords_rectangle(0,
                                           tg.height_dielectric +
                                           tg.height_copper + tg.height_gap,
                                           tg.radius_pcb,
                                           tg.height_pcb_core +
                                           2 * tg.height_pcb_prepreg +
                                           2 * tg.height_copper)
    else:
        print('Unknown number of layers on primary side')
        femm.closefemm()
        return (0, 0)

    # coordinates of the PCB on the secondary side
    if tg.layers_secondary == 1 or tg.layers_secondary == 2:
        coords_pcb_sec = coords_rectangle(0,
                                          -tg.height_copper - tg.height_gap,
                                          tg.radius_pcb,
                                          -tg.height_pcb_core)
    elif tg.layers_secondary == 4:
        coords_pcb_sec = coords_rectangle(0,
                                          -tg.height_copper - tg.height_gap,
                                          tg.radius_pcb,
                                          -tg.height_pcb_core -
                                          2 * tg.height_pcb_prepreg -
                                          2 * tg.height_copper)
    else:
        print('Unknown number of layers on secondary side')
        femm.closefemm()
        return (0)
    femm.ei_drawrectangle(*coords_pcb_prim)
    femm.ei_drawrectangle(*coords_pcb_sec)
    label_dict['PCB_prim'] = coords_pcb_prim
    label_dict['PCB_sec'] = coords_pcb_sec


def add_isolation(tg):
    '''Add isolation disc and surrounding gel or liquid.'''
    # coordinates of isolation disc
    coords_isolation = coords_rectangle(0, 0, tg.radius_dielectric,
                                        tg.height_dielectric)
    # coordinates of gel or liquid surrounding the transformer
    coords_gel = coords_rectangle(0, -tg.height_gel, tg.radius_gel,
                                  2 * tg.height_gel + tg.height_dielectric)
    femm.ei_drawrectangle(*coords_isolation)
    femm.ei_drawrectangle(*coords_gel)


def add_magnetics(tg, label_dict):
    ''' '''
    # TODO: add function
    for m in tg.magnetics:
        # draw boxes
        if not m.polarity:  # primary side
            start_z = label_dict['prim' + str(tg.layers_primary*tg.turns_primary-1)][1] + m.distance # take copper into accoutns
            coords = coords_rectangle(m.radius_inner, start_z, m.radius_outer - m.radius_inner, m.thickness)
        else:   # secondary side
            start_z = label_dict['sec' + str(tg.layers_secondary*tg.turns_secondary-1)][1] - m.distance # take copper into accoutns
            coords = coords_rectangle(m.radius_inner, start_z, m.radius_outer - m.radius_inner, -m.thickness)
        femm.ei_drawrectangle(*coords)

        # add label
        label_coord = (np.average((coords[0], coords[2])),
                       np.average((coords[1], coords[3])))
        femm.ei_addblocklabel(*label_coord)
        femm.ei_selectlabel(*label_coord)
        femm.ei_setblockprop(m.material, 1, 0, 'None')
        femm.ei_clearselected()
        # TODO add label to label dictionnary
        # set boundary conditions
        if m.polarity:
            set_conductor_boundry(coords, 0)  # secondary
        else:
            set_conductor_boundry(coords, 1)  # primary


def add_block_labels(tg, label_dict):
    '''Add block labels to pcbs, air, isolation disc, isolation gel/liquid.'''
    # Add block labels
    # ei seteditmode(editmode)
    # From manual: Sets the current editmode to:
    # – "nodes" - nodes
    # – "segments" - line segments
    # – "arcsegments" - arc segments
    # – "blocks" - block labels
    # – "group" - selected group
    femm.ei_seteditmode('blocks')

    # ei addblocklabel(x,y)
    # From manual: Add a new block label at (x,y)

    # labels for the PCBs
    coords = [2, tg.height_dielectric + tg.height_pcb_core / 2.]
    femm.ei_addblocklabel(*coords)
    label_dict['pcb_prim_label'] = coords

    coords = [2, - tg.height_pcb_core / 2.]
    femm.ei_addblocklabel(*coords)
    label_dict['pcb_sec_label'] = coords

    # label for the surrounding air
    coords = [2, tg.height_dielectric * 2 + tg.height_gel]
    femm.ei_addblocklabel(*coords)
    label_dict['air'] = coords

    # label for the dilectric gel or liquid surrounding the transformer
    coords = [tg.radius_pcb + 2, tg.height_dielectric + tg.height_gel / 2.]
    femm.ei_addblocklabel(*coords)
    label_dict['gel'] = coords

    # label for the isolation disc
    coords = [2, tg.height_dielectric / 2.]
    femm.ei_addblocklabel(*coords)
    label_dict['insulation'] = coords

    # Set material type for all blocks
    # ei setblockprop("blockname", automesh, meshsize, group) Set the selected
    # block labels to have the properties: Block property "blockname".
    # automesh: 0 = mesher defers to mesh size constraint defined in meshsize,
    # 1 = mesher automatically chooses the mesh density. meshsize: size
    # constraint on the mesh in the block marked by this label. A member of
    # group number group

    femm.ei_selectlabel(*label_dict['pcb_prim_label'])
    femm.ei_selectlabel(*label_dict['pcb_sec_label'])
    femm.ei_setblockprop('fr4', 1, 0, 'None')
    femm.ei_clearselected()

    femm.ei_selectlabel(*label_dict['air'])
    femm.ei_setblockprop('air', 1, 0, 'None')
    femm.ei_clearselected()

    femm.ei_selectlabel(*label_dict['gel'])
    femm.ei_setblockprop('midel', 1, 0, 'None')
    femm.ei_clearselected()

    femm.ei_selectlabel(*label_dict['insulation'])
    femm.ei_setblockprop(tg.material_dielectric.lower(), 1, 0, 'None')
    femm.ei_clearselected()


def modify_tracks(tg, label_dict, segment_angle=1, segment_length=0.004):
    '''Modify tracks according to the FancyTrack class. This includes rounding
    of corners and elongation of individual tracks.

    Args:
        tg (:obj:'TransformerGeometry'): contains all geometry information
        needed to build conductors in the problem.
        label_dict (dict): stores label coordinates
    '''
    for track in tg.tracks:
        if track.elongation != 0:
            _modify_tracks_elongate(tg, track, label_dict)
        if track.rounding:
            coords, vectors, labels = _get_corner_coord(track, tg, label_dict)
            _modify_tracks_round(track, coords, tg.radius_corner,
                                 segment_angle)

    _increase_pcb_mesh(tg, 3*tg.height_copper, segment_length, label_dict)


def _modify_tracks_round(track, coords, radius, segment_angle):
    # Define parameters of the rounding
    if track.polarity == 0:
        volt_potential = 'high'
    else:
        volt_potential = 'zero'
    for c in coords:
        femm.ei_createradius(c[0], c[1], radius)
        femm.ei_selectarcsegment(*c)
        femm.ei_setarcsegmentprop(segment_angle, '<None>', 0, 0,
                                  volt_potential)
        femm.ei_clearselected()


def _modify_tracks_elongate(tg, track, label_dict):
    ''' '''
    if not track.polarity:
        label = 'prim'
        turns_tot = tg.turns_primary
        cu_width = tg.width_copper[0]
    else:
        label = 'sec'
        turns_tot = tg.turns_secondary
        cu_width = tg.width_copper[1]
    label = label + str(track.layer * turns_tot + track.turn)
    coords = label_dict[label]

    # elongate inner and/or outer line segment.
    if track.side_h in ['inner', 'both']:
        coords_inner = coords + np.array([-cu_width / 2, 0])
        femm.ei_selectsegment(*coords_inner)
        femm.ei_movetranslate2(-track.elongation, 0, 1)
        femm.ei_clearselected()

    if track.side_h in ['outer', 'both']:
        coords_outer = coords + np.array([cu_width / 2, 0])
        femm.ei_selectsegment(*coords_outer)
        femm.ei_movetranslate2(track.elongation, 0, 1)
        femm.ei_clearselected()


def _increase_pcb_mesh(tg, cut_in, segment_length, label_dict):
    '''Decrease segment size to increase mesh along edge of the innermost PCB
    edge where the field is high.

    Args:
        tg (:obj:'TransformerGeometry'): contains all geometry information
        cut_in (float): how far along the pcb edge from the conductor where the
        segment length is set
        dz (float): copper height of conductors
        segment_length (float): specifiy segment length in femm
        label_dict (dict): stores label coordinates
    '''
    # Define parameters of the rounding
    cu_height = tg.height_copper
    cu_width = tg.width_copper
    segment_length = 0.004
    pcb_segment = []
    pcb_node = []

    # Find primary side corners and add new labels for the rounded edge
    edge0 = (np.array(label_dict['prim0']) +
             np.array((-cu_width[0],
                       cu_height)) / 2.0)
    edge1 = (np.array(label_dict['prim' + str(tg.turns_primary - 1)]) +
             np.array((cu_width[0],
                       cu_height)) / 2.0)
    # Find secondary side corners and add new labels for the rounded edge
    edge2 = (np.array(label_dict['sec0']) +
             np.array((-cu_width[1],
                       -cu_height)) / 2.0)
    edge3 = (np.array(label_dict['sec' + str(tg.turns_secondary - 1)]) +
             np.array((cu_width[1],
                       -cu_height)) / 2.0)

    pcb_node.append(edge0 + np.array((-cut_in, 0)))
    pcb_segment.append(pcb_node[-1] + np.array((cut_in / 2, 0)))
    pcb_node.append(edge1 + np.array((cut_in, 0)))
    pcb_segment.append(pcb_node[-1] + np.array((-cut_in / 2, 0)))
    pcb_node.append(edge2 + np.array((-cut_in, 0)))
    pcb_segment.append(pcb_node[-1] + np.array((cut_in / 2, 0)))
    pcb_node.append(edge3 + np.array((cut_in, 0)))
    pcb_segment.append(pcb_node[-1] + np.array((-cut_in / 2, 0)))

    for idx, node in enumerate(pcb_node):
        femm.ei_addnode(*node)

    for segment in pcb_segment:
        femm.ei_selectsegment(*segment)
        femm.ei_setsegmentprop('<None>', segment_length, 0, 0, 0, '<None>')
        femm.ei_clearselected()


def add_guard(geometry, guard, label_dict):
    ''' Add guard ring(s) in the problem

    Args:
        tg (TransformerGeometry): contains all geometry information needed
        to build conductors in the problem.
        guard (:obj:'list: of :obj:'Guard'): contains guard geometry
        label_dict (dict): stores label coordinates
    '''

    if not isinstance(guard, (list, tuple)):
        guard = [guard, ]
    for idx, g in enumerate(guard):
        if g.gtype == 'Ring':
            draw_guard_ring(g.distance, g.gap, g.radius, g.polarity,
                            'guard' + str(idx), label_dict)
        elif g.gtype == 'Trench':
            draw_guard_trench(geometry, g, 'trench' + str(idx), label_dict)


def draw_guard_ring(distance, gap, radius, polarity, label_name, label_dict):
    '''Draw circular shape of guard ring in the problem and add a label.

    Args:
        distance (float): distance in r-plane from edge of PCB to edge of
        guard ring
        gap (float): distance in z-plane from edge of PCB to edge of guard ring
        radius (float): radius of the guard ring conductor
        polarity (int): attach the gaurd ring to high side or low side.
        Polarity higher than zero sets the guard ring to the upper side of the
        dielectric. Else the guard ring is set to lower side.
        label_name (str): name of this conductor that is stored in the label
        dictionary
        label_dict (dict): stores label coordinates
    '''

    max_angle = 1  # arc_section angle, affects finess of mesh.
    # high or low side?
    if polarity > 0:
        origin = label_dict['PCB_prim']
        origin = np.array((origin[2], origin[1]))
        point1 = origin + np.array((radius + distance, gap))
        point2 = point1 + np.array((0, 2 * radius))
        in_conductor = 'high'
    else:
        origin = label_dict['PCB_sec']
        origin = np.array((origin[2], origin[1]))
        point1 = origin + np.array((radius + distance, -gap))
        point2 = point1 + np.array((0, -2 * radius))
        in_conductor = 'zero'
    label_dict[label_name + 'surface'] = point1

    # ei_drawarc(x1,y1,x2,y2,angle,maxseg)
    # From manual: Adds nodes at (x1,y1) and (x2,y2) and adds an arc of the
    # specified angle and discretization connecting the nodes.
    femm.ei_drawarc(point1[0], point1[1], point2[0], point2[1], 180, max_angle)
    femm.ei_addarc(point2[0], point2[1], point1[0], point1[1], 180, max_angle)

    # Add label and block property to guard ring
    label_coord = np.average((point1, point2), axis=0)
    femm.ei_addblocklabel(*label_coord)
    femm.ei_selectlabel(*label_coord)
    femm.ei_setblockprop('<No Mesh>', 1, 0, 'guards')
    femm.ei_clearselected()
    label_dict[label_name] = label_coord

    # Set boundry based on coordinates
    # ei_setarcsegmentprop(maxsegdeg, ’propname’, hide, group, ’inconductor’)
    # From manual: Set the selected arc segments to:
    # Meshed with elements that span at most maxsegdeg degrees per element
    # Boundary property ’propname’
    # hide: 0 = not hidden in post-processor, 1 == hidden in post processor
    # A member of group number group
    # A member of the conductor specified by the string ’inconductor’. If the
    # segment is not part of a conductor, this parameter can be specified
    # as ’<None>’.
    femm.ei_selectarcsegment(label_coord[0] - radius, label_coord[1])
    femm.ei_selectarcsegment(label_coord[0] + radius, label_coord[1])
    femm.ei_setarcsegmentprop(max_angle, 'None', 0, 'None', in_conductor)
    femm.ei_clearselected()


def draw_guard_trench(tg, trench, label_name, label_dict):
    '''Draw a trench on the outside of turn closest to the dielectric.

    Args:
        tg (:obj:'TransformerGeometry'): contains all geometry information
        trench (:obj:'GuardTrench'): trench guard geometry
        label_name (str): label that identifies the guard in the dict
        label_dict (dict): stores label coordinates
    '''

    conductor_height = tg.height_copper
    conductor_width = tg.width_copper
    max_angle = 2
    elementsize = 0.05
    if trench.inner:
        sign = -1
        turns = 1
    else:
        sign = 1
        turns = tg.turns_primary

    # find starting point
    if trench.polarity:
        node0 = (np.array(label_dict['prim' + str(turns - 1)]) +
                 np.array((conductor_width[0]*sign, conductor_height))/2.0 +
                 np.array((trench.distance * sign, 0)))
        in_conductor = 'high'
        pm = 1  # polarity multiplier
    else:
        node0 = (np.array(label_dict['sec' + str(turns - 1)]) +
                 np.array((conductor_width[1]*sign, -conductor_height))/2.0 +
                 np.array((trench.distance * sign, 0)))
        in_conductor = 'zero'
        pm = -1  # polarity multiplier
    # draw four/five points and lines between points 1 to last to 0
    nodes = [(0, 0), (trench.width*sign, 0),
             (trench.width*sign, trench.depth * pm), (0, trench.depth * pm)]
    if trench.drill_angle > 90:
        alpha = (1.0 - trench.drill_angle/180.) * pi
        depth_added = trench.width / (2.0 * np.tan(alpha))
        nodes.insert(-1, (trench.width/2.0*sign,
                          (trench.depth + depth_added) * pm))
    nodes = [node0 + np.array(n) for n in nodes]
    for idx, n in enumerate(nodes):
        femm.ei_addnode(*n)
        label_dict[label_name + 'node'+str(idx)] = n
    ind = range(1, len(nodes)) + [0, ]
    segment_gen = (np.concatenate((nodes[ind[idx]], nodes[ind[idx+1]]))
                   for idx in range(len(nodes)-1))
    for n in segment_gen:
        femm.ei_addsegment(*n)
        femm.ei_selectsegment(*(np.mean((n[:2], n[2:]), axis=0)))
    femm.ei_setsegmentprop('None', elementsize, 0, 0, 'None', in_conductor)
    femm.ei_clearselected()
    # remove existing segment of the PCB
    femm.ei_selectsegment(*(np.mean((nodes[0:2]), axis=0)))
    femm.ei_deleteselectedsegments()
    # draw circular top
    if (pm == 1 and sign == 1) or (pm == -1 and sign == -1):
        femm.ei_addarc(nodes[0][0], nodes[0][1], nodes[1][0], nodes[1][1],
                       trench.outer_angle, max_angle)
    else:
        femm.ei_addarc(nodes[1][0], nodes[1][1], nodes[0][0], nodes[0][1],
                       trench.outer_angle, max_angle)

    # set conductor properties on the boundries
    femm.ei_selectarcsegment(*(np.mean((nodes[0:2]), axis=0)))
    femm.ei_setarcsegmentprop(max_angle, 'None', 0, 'None', in_conductor)

    # add label and set material
    label_dict[label_name] = (nodes[0] +
                              np.array((trench.width/2.*sign,
                                        trench.depth*pm/2.)))
    femm.ei_addblocklabel(*label_dict[label_name])
    femm.ei_selectlabel(*label_dict[label_name])
    femm.ei_setblockprop(trench.material, 1, 0, 'None')
    femm.ei_clearselected()


def draw_field_contour(coords_start, coords_end, filename=None):
    '''Make (and save) countour plot.'''
    # eo_seteditmode(mode)
    # From manual: Sets themode of the postprocessor to point, contour, or area
    # mode. Valid entries for mode are "point", "contour", and "area".
    femm.eo_seteditmode('contour')

    # eo_selectpoint(x,y)
    # From manual: Adds a contour point at the closest input point to (x,y).
    # If the selected point and a previously selected point lie at the ends of
    # an arcsegment, a contour is added that traces along the arcsegment. The
    # selectpoint command has the same functionality as the left-button-click
    # contour point selection when the program is running in interactive mode.
    femm.eo_selectpoint(*coords_start)
    femm.eo_selectpoint(*coords_end)

    # eo_makeplot(PlotType,NumPoints,Filename,FileFormat)
    # From manual: Allows Octave to access to the X-Y plot routines. If only
    # PlotType or only PlotType and NumPoints are specified, the command is
    # interpreted as a request to plot the requested plot type to the screen.
    # If, in addition, the Filename parameter is specified, the plot is instead
    # written to disk to the specified file name as an extended metafile. If
    # the FileFormat parameter is also, the command is instead interpreted as a
    # command to write the data to disk to the specfied file name, rather than
    # display it to make a graphical plot. Valid entries for PlotType are:
    # PlotType Definition
    # 0 V (Voltage)
    # 1 |D| (Magnitude of flux density)
    # 2 D . n (Normal flux density)
    # 3 D . t (Tangential flux density)
    # 4 |E| (Magnitude of field intensity)
    # 5 E . n (Normal field intensity)
    # 6 E . t (Tangential field intensity)
    # Valid file formats are:
    # FileFormat Definition
    # 0 Multi-column text with legend
    # 1 Multi-column text with no legend
    # 2 Mathematica-style formatting

    # eo_clearcontour
    # From manual: Clear a prevously defined contour
    if filename:
        femm.eo_makeplot(4, 5000, filename, 0)
        femm.eo_clearcontour()
    else:
        femm.eo_makeplot(4, 5000)


def get_peak_field(tg, label_dict, segment_angle=1, distance=None):
    ''' Get value of electric field for all rounded corners that are given as
    fancy tracks of the geometry. The field value is ziped toghether with a
    text label that indicates which corner.'''
    e_field = []

    for t in tg.tracks:
        if t.rounding:
            coords, vectors, labels = _get_corner_coord(t, tg, label_dict)
            for c, v, l in zip(coords, vectors, labels):
                e_field.append((l, _get_peak_field(c, v, tg.radius_corner,
                                                   segment_angle, distance)))
    return e_field


def _get_peak_field(coord, vector, radius, segment_angle, distance):
    ''' Get the E-field value a given distance away from a rounded corner. If
    no distance is given, the distance is set to be 5 times the mesh size of
    the arc segment that makes up the corner.'''
    if distance is None:
        distance = 5*pi*radius*segment_angle/180
    coord_shift = np.array([-np.sqrt(2)/2 * (radius + distance)]*2) + radius
    coord_point = coord - vector * coord_shift
    x = femm.eo_getpointvalues(*coord_point)
    return np.sqrt(x[3]**2 + x[4]**2)


def _get_corner_coord(track, tg, label_dict):
    coord_corner = []
    vector_corner = []
    label_corner = []

    if track.rounding is not None:
        if not track.polarity:
            label = 'prim'
            turns_tot = tg.turns_primary
            cu_width = tg.width_copper[0]
        else:
            label = 'sec'
            turns_tot = tg.turns_secondary
            cu_width = tg.width_copper[1]
        label = label + str(track.layer * turns_tot + track.turn)
        coord = label_dict[label]
        cu_height = tg.height_copper
        elongation = track.elongation

        # side
        if track.side_h in ['inner', 'both']:
            vector_h = -1
            label_base = label + '_inner'
            if track.side_v in ['towards', 'both']:
                vector_v = -1 + 2 * track.polarity
                vector = np.array([vector_h, vector_v])

                x = np.array((cu_width / 2.0 * vector[0] - elongation,
                              cu_height / 2.0 * vector[1]))
                coord_corner.append(coord + x)
                vector_corner.append(vector)
                label_corner.append(label_base + '_towards')
            if track.side_v in ['away', 'both']:
                vector_v = 1 - 2 * track.polarity
                vector = np.array([vector_h, vector_v])

                x = np.array((cu_width / 2.0 * vector[0] - elongation,
                              cu_height / 2.0 * vector[1]))
                coord_corner.append(coord + x)
                vector_corner.append(vector)
                label_corner.append(label_base + '_away')
        if track.side_h in ['outer', 'both']:
            label_base = label + '_outer'
            vector_h = 1
            if track.side_v in ['towards', 'both']:
                vector_v = -1 + 2 * track.polarity
                vector = np.array([vector_h, vector_v])

                x = np.array((cu_width / 2.0 * vector[0] + elongation,
                              cu_height / 2.0 * vector[1]))
                coord_corner.append(coord + x)
                vector_corner.append(vector)
                label_corner.append(label_base + '_towards')
            if track.side_v in ['away', 'both']:
                vector_v = 1 - 2 * track.polarity
                vector = np.array([vector_h, vector_v])

                x = np.array((cu_width / 2.0 * vector[0] + elongation,
                              cu_height / 2.0 * vector[1]))
                coord_corner.append(coord + x)
                vector_corner.append(vector)
                label_corner.append(label_base + '_away')
    return [coord_corner, vector_corner, label_corner]


def _skip_last(iterator):
    '''Generator that skips the last entry. Used when reading output from
    contour plots as the final line of the output file is empty.'''
    prev = next(iterator)
    for item in iterator:
        yield prev
        prev = item


def read_field_contour_file(filename):
    '''Read csv file of a contour plot. Returns the text string for the units
    in the x and y axis, and the data as two lists [x] [y].'''
    with open(filename, 'r') as f:
        unit_xaxis = f.readline()
        unit_yaxis = f.readline()
        data = []
        for row in _skip_last(f):
            a, b, c = row.split('\t')
            data.append([float(a), float(b)])
        return unit_xaxis, unit_yaxis, list(zip(*data))


def get_field_contour_max(filename):
    '''Read contour out file and return the maximum value.'''
    unitx, unity, values = read_field_contour_file(filename)
    return np.max(values[1])


def save_field_contour(conductor_height, label_dict, guard=None):
    '''Make field contour plots of critical points in the transformer geometry;
    edges and guards.'''
    edge_radius = conductor_height * .2
    # edge0 ++ == upper inner conductor
    # edge1 -+ == upper outer conductor
    # edge2 +- == lower inner conductor
    # edge3 -- == lower outer conductor
    sign = np.array(((1, 1), (-1, 1), (1, -1), (-1, -1)))
    femm.eo_seteditmode('contour')

    for idx in range(4):
        coords = label_dict['edge' + str(idx)]
        contour_start = coords
        contour_stop = coords + (np.array((edge_radius, edge_radius)) *
                                 sign[idx])
        femm.eo_addcontour(*contour_start)
        femm.eo_addcontour(*contour_stop)
        femm.eo_makeplot(4, 5000, '\\output_field_edge' + str(idx) + '.csv', 0)
        femm.eo_clearcontour()

    if guard:
        for idx in range(int(len(guard) / 2)):
            contour_top = np.array(label_dict['guard' + str(idx*2)])
            contour_top += np.array((0, guard[idx*2][2]))
            contour_bottom = np.array(label_dict['guard' + str(idx*2+1)])
            contour_bottom += np.array((0, -guard[idx*2+1][3]))
            fname = '\\output_field_guard' + str(idx) + '.csv'
            draw_field_contour(contour_top, contour_bottom, filename=fname)


def plot_field_contour(filename):
    ''' plot contour output file data '''
    f = plt.figure()
    for f in list(filename):
        unitx, unity, values = read_field_contour_file(f)
        plt.plot(values[0], values[1], linewidth=2)
        plt.xlabel(unitx[unitx.find(':') + 2:])
        plt.ylabel(unity[unity.find(':') + 2:])
    return f


def set_view(view, guard, label_dict):
    '''Change view (and optionally save files).'''
    for v, l in get_zoom_coords_corners(view, label_dict):
        femm.eo_zoom(*v)
        if view.filename:
            femm.eo_savebitmap(view.filedir + view.filename + l + '.bmp')
    if guard:
        for v, l in get_zoom_coords_guard(view, guard, label_dict):
            femm.eo_zoom(*v)
            if view.filename:
                femm.eo_savebitmap(view.filedir + view.filename + l + '.bmp')
    else:
        # end the view on upper side outer corner if there is no guard,
        # as this is usually the most critical part of the design
        femm.eo_zoom(*get_zoom_coords_corners(view, label_dict)[1][0])


def get_zoom_coords_corners(view, label_dict):
    '''Set view for edges.'''
    span = np.array((view.span_r, view.span_z)) / 2.0
    zoom_coords = []
    labels = []
    for label in ['edge'+str(idx) for idx in range(4)]:
        center_point = np.array(label_dict[label])
        labels.append(label)
        zoom_coords.append(np.concatenate((center_point - span,
                                           center_point + span)))
    return zip(zoom_coords, labels)


def get_zoom_coords_guard(view, guard, label_dict):
    '''Set view for guards.'''
    zoom_coords = []
    labels = []
    if not isinstance(guard, (list, tuple)):
        guard = [guard, ]
    for idx, g in enumerate(guard):
        span = np.array((view.span_r, view.span_z))
        if g.gtype == 'Ring':
            label = 'guard{}surface'
            if (span < np.array((2*g.radius, g.radius))).all():
                span = np.array((2*g.radius, g.radius)) / 2.0
        elif g.gtype == 'Trench':
            label = 'trench{}'
            if (span < 2 * np.array((g.width, g.depth))).all():
                span = np.array((g.width, g.width*9./16.))
        else:
            span = span / 2.0
        labels.append(label.format(idx))
        center_point = np.array(label_dict[label.format(idx)])
        zoom_coords.append(np.concatenate((center_point - span,
                                           center_point + span)))
    return zip(zoom_coords, labels)


def calc_field_distribution(tg, voltage_high=1, view=None, **kwargs):
    ''' Setup of electro-static problem in femm to calculate capacitance
    between primary and secondary side of planar transformer. Additionally,
    electric field can be investigated

    Args:
        tg (:obj:'TransformerGeometry'): tg contains all geometry information
        needed to solve the problem.
        voltage_high (float): Voltage potential of the primary side.
        view (:obj:'View'): Save images of edges and guards.

    Returns:
        capacitance (float): calculated capacitance between primary and
        secondary side of planar transformer.
    '''

    etiquettes_dict = {}  # Dictionary to store coordinates of nodes

    # initialitiation of the electro-static problem
    boundary_radius = 2 * tg.radius_dielectric
    initial_setup(boundary_radius, voltage_high, **kwargs)
    if 'segment_angle' in kwargs:
        segment_angle = kwargs('segment_angle')
    else:
        segment_angle = 1

    # draw conductors on primary and secondary sides
    add_conductors(tg, etiquettes_dict)
    add_pcbs(tg, etiquettes_dict)
    add_isolation(tg)
    add_magnetics(tg, etiquettes_dict)
    add_block_labels(tg, etiquettes_dict)

    modify_tracks(tg, etiquettes_dict, segment_angle)

    if tg.guard:
        add_guard(tg, tg.guard, etiquettes_dict)

    # mi zoomnatural()
    # From manual: zooms to a “natural” view with sensible extents.
    femm.ei_zoomnatural()

    # Save geometry file
    femm.ei_saveas('field_calculations.fee')

    # Meshing and analysis
    # From manual: Note that it is not necessary to run mesh before performing
    # an analysis, as mi_analyze() will make sure the mesh is up to date before
    # running an analysis.

    # mi analyze(flag)
    # From manual: runs fkern to solve the problem. The flag parameter controls
    # whether the fkern window is visible or minimized. For a visible window,
    # either specify no value for flag or specify 0. For a minimized window,
    # flag should be set to 1.
    femm.ei_analyze(1)

    # Post-processor
    femm.ei_loadsolution()

    # eo_getconductorproperties("conductor")
    # From manual: Properties are returned for the conductor property named
    # ”conductor”. Two values are returned: The voltage of the specified
    # conductor, and the charge carried on the specified conductor.

    # Calculate the capacitance
    circuit_properties = femm.eo_getconductorproperties('zero')
    charge = -circuit_properties[1]  # get charge on secondary side conductors
    capacitance = charge / voltage_high

    # eo_showdensityplot(legend,gscale,type,upper D,lower D)
    # From manual: Shows the flux density plot with options:
    # legend Set to 0 to hide the plot legend or 1 to show the plot legend.
    # gscale Set to 0 for a colour density plot or 1 for a grey scale density
    # plot.
    # type Sets the type of density plot. A value of 0 plots voltage, 1 plots
    # the magnitude of D, and 2 plots the magnitude of E
    # upper D Sets the upper display limit for the density plot.
    # lower D Sets the lower display limit for the density plot.
    if 'plot_field_max' in kwargs:
        femm.eo_refreshview()  # need to refresh before a change can be made
        femm.eo_showdensityplot(1, 0, 2, kwargs['plot_field_max'] * 20 / 19, 0)

    # Save field lines along outer point of the PCB conductors.
    if view:
        set_view(view, tg.guard, etiquettes_dict)

    # Find high field values of corners
    e_field_corners = get_peak_field(tg, etiquettes_dict)

    if kwargs.get('close') is True:
        femm.closefemm()

    return capacitance, e_field_corners
