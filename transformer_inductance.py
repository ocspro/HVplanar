#!/usr/bin/env python
# -*- coding: utf-8 -*-

from numpy import average

import femm


def initial_setup(boundary, currents, **kwargs):
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
    femm.newdocument(0)

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
    if 'frequency' in kwargs:
        freq = kwargs['frequency']
    else:
        freq = 6.78e6
    if 'precision' in kwargs:
        precision = kwargs['precision']
    else:
        precision = 5e-9
    if 'min_angle' in kwargs:
        min_angle = kwargs['min_angle']
    else:
        min_angle = 30
    femm.mi_probdef(freq, 'millimeters', 'axi', precision, 50, min_angle)

    # Circuit parameters
    # mi addcircprop("circuitname", i, circuittype)
    # From manual: adds a new circuit property with name "circuitname" with a
    # prescribed current, i. The circuittype parameter is 0 for a
    # parallel-connected circuit and 1 for a series-connected circuit.

    # The currents in the primary and secondary circuits are set
    I1, I2 = currents
    femm.mi_addcircprop('phase_prim', I1, 1)
    femm.mi_addcircprop('phase_sec', I2, 1)

    # Add materials properties used in the simulation
    # mi addmaterial("materialname", mu x, mu y, H c, J, Cduct, Lam d,
    # Phi hmax, lam fill, LamType, Phi hx, Phi hy,NStrands,WireD)
    # From manual: adds a newmaterial with called "materialname" with the
    # material properties:
    # – mu x Relative permeability in the x- or r-direction.
    # – mu y Relative permeability in the y- or z-direction.
    # – H c Permanent magnet coercivity in Amps/Meter.
    # – J Real Applied source current density in Amps/mm2.
    # – Cduct Electrical conductivity of the material in MS/m.
    # – Lam d Lamination thickness in millimeters.
    # – Phi hmax Hysteresis lag angle in degrees, used for nonlinear BH curves.
    # – Lam fill Fraction of the volume occupied per lamination that is
    # actually filled with iron (Note that this parameter defaults to 1 the
    # femm preprocessor dialog box because, by default, iron completely fills
    # the volume)
    # – Lamtype Set to
    # ? 0 – Not laminated or laminated in plane
    # ? 1 – laminated x or r
    # ? 2 – laminated y or z
    # ? 3 – Magnet wire
    # ? 4 – Plain stranded wire
    # ? 5 – Litz wire
    # ? 6 – Square wire
    # – Phi hx Hysteresis lag in degrees in the x-direction for linear problems
    # – Phi hy Hysteresis lag in degrees in the y-direction for linear problems
    # – NStrands Number of strands in the wire build. Should be 1 for Magnet or
    # Square wire.
    # – WireD Diameter of each wire constituent strand in millimeters.
    femm.mi_addmaterial('air', 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    femm.mi_addmaterial('fr4', 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    femm.mi_addmaterial('copper', 1, 1, 0, 0, 58, 0, 0, 0, 0, 0, 0)
    femm.mi_addmaterial('polysterimide', 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    femm.mi_addmaterial('teflon', 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    femm.mi_addmaterial('silgel', 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    if 'material' in kwargs:
        for m in kwargs['material']:
            femm.mi_addmaterial(*m)

    # Boundary condition
    # ei makeABC(n,R,x,y,bc)
    # From manual: creates a series of circular shells that emulate the
    # impedance of an unbounded domain (i.e. an Improvised Asymptotic Boundary
    # Condition). The n parameter contains the number of shells to be used
    # (should be between 1 and 10), R is the radius of the solution domain, and
    # (x,y) denotes the center of the solution domain. The bc parameter should
    # be specified as 0 for a Dirichlet outer edge or 1 for a Neumann outer
    # edge. If the function is called without all the parameters, the function
    # makes up reasonable values for the missing parameters.
    femm.mi_makeABC(7, boundary, 0, 0, 0)


def coords_of_rectangle(x0, y0, dx, dy):
    ''' Return array with the coordinates of the two opposite corners of a
    rectangle.'''
    return (x0, y0, x0+dx, y0+dy)


def get_z_coord_copper_layer(g):
    ''' Calculate z coordinate of copper layers depending on the number of
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


def draw_conductor(coords, in_conductor, label_name, label_dict):
    '''Draw a square conductor, set a block label in the middle of the
    conductor and add label properties.'''
    femm.mi_drawrectangle(*coords)

    label_coord = (average((coords[0], coords[2])),
                   average((coords[1], coords[3])))
    femm.mi_addblocklabel(*label_coord)
    femm.mi_selectlabel(*label_coord)
    if in_conductor == 1:
        femm.mi_setblockprop('copper', 1,  0, 'phase_prim', 0, 1, 1)
    elif in_conductor == 0:
        femm.mi_setblockprop('copper', 1,  0, 'phase_sec', 0, 1, 1)
    femm.mi_clearselected()
    label_dict[label_name] = label_coord


def add_conductors(g, label_dict):
    '''Add all conductors of primary and secondary side pcbs

    Args:
        g (:obj:'TransformerGeometry'): g contains all geometry information
        needed to build conductors in the problem.
    '''

    z_copper = get_z_coord_copper_layer(g)
    # adding conductors on the primary side
    for c in range(g.layers_primary):
        for i in range(g.turns_primary):
            coords_conductor = coords_of_rectangle(g.radius_inner_track[0] +
                                                   (g.width_between_tracks[0] +
                                                    g.width_copper[0]) * i,
                                                   z_copper[0][c],
                                                   g.width_copper[0],
                                                   g.height_copper)
            draw_conductor(coords_conductor,
                           1,
                           'prim' + str(i + c * g.turns_primary),
                           label_dict)
    # adding conductors on the secondary side
    for c in range(g.layers_secondary):
        for i in range(g.turns_secondary):
            coords_conductor = coords_of_rectangle(g.radius_inner_track[1] +
                                                   (g.width_between_tracks[1] +
                                                    g.width_copper[1]) * i,
                                                   z_copper[1][c],
                                                   g.width_copper[1],
                                                   -g.height_copper)
            draw_conductor(coords_conductor,
                           0,
                           'sec' + str(i + c * g.turns_secondary),
                           label_dict)


def add_pcbs(tg):
    '''Add primary and secondary side pcbs.'''
    # coordinates of the PCB on the primary side
    if tg.layers_primary == 1 or tg.layers_primary == 2:
        coords_pcb_prim = coords_of_rectangle(0,
                                              tg.height_dielectric +
                                              tg.height_copper + tg.height_gap,
                                              tg.radius_pcb,
                                              tg.height_pcb_core)
    elif tg.layers_primary == 4:
        coords_pcb_prim = coords_of_rectangle(0,
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
        coords_pcb_sec = coords_of_rectangle(0,
                                             -tg.height_copper - tg.height_gap,
                                             tg.radius_pcb,
                                             -tg.height_pcb_core)
    elif tg.layers_secondary == 4:
        coords_pcb_sec = coords_of_rectangle(0,
                                             -tg.height_copper - tg.height_gap,
                                             tg.radius_pcb,
                                             -tg.height_pcb_core -
                                             2 * tg.height_pcb_prepreg -
                                             2 * tg.height_copper)
    else:
        print('Unknown number of layers on secondary side')
        femm.closefemm()
        return (0, 0)
    femm.mi_drawrectangle(*coords_pcb_prim)
    femm.mi_drawrectangle(*coords_pcb_sec)


def add_isolation(tg):
    # coordinates of isolation disc
    coords_isolation = coords_of_rectangle(0, 0, tg.radius_dielectric,
                                           tg.height_dielectric)
    # coordinates of gel or liquid surrounding the transformer
    coords_gel = coords_of_rectangle(0, -tg.height_gel, tg.radius_gel,
                                     2 * tg.height_gel + tg.height_dielectric)
    femm.mi_drawrectangle(*coords_isolation)
    femm.mi_drawrectangle(*coords_gel)


def add_block_labels(tg, etiquettes_dict):
    '''Add block labels to pcbs, air, isolation disc, isolation gel/liquid.'''
    # Add block labels
    # ei seteditmode(editmode)
    # From manual: Sets the current editmode to:
    # – "nodes" - nodes
    # – "segments" - line segments
    # – "arcsegments" - arc segments
    # – "blocks" - block labels
    # – "group" - selected group
    femm.mi_seteditmode('blocks')

    # ei addblocklabel(x,y)
    # From manual: Add a new block label at (x,y)

    # labels for the PCBs
    coords = [2, tg.height_dielectric + tg.height_pcb_core / 2.]
    femm.mi_addblocklabel(*coords)
    etiquettes_dict['pcb_prim'] = coords

    coords = [2, - tg.height_pcb_core / 2.]
    femm.mi_addblocklabel(*coords)
    etiquettes_dict['pcb_sec'] = coords

    # label for the surrounding air
    coords = [2, tg.height_dielectric * 2 + tg.height_gel]
    femm.mi_addblocklabel(*coords)
    etiquettes_dict['air'] = coords

    # label for the dilectric gel or liquid surrounding the transformer
    coords = [tg.radius_pcb + 2, tg.height_dielectric + tg.height_gel / 2.]
    femm.mi_addblocklabel(*coords)
    etiquettes_dict['gel'] = coords

    # label for the isolation disc
    coords = [2, tg.height_dielectric / 2.]
    femm.mi_addblocklabel(*coords)
    etiquettes_dict['isolant'] = coords

    # Set material type for all blocks
    # mi_setblockprop("blockname", automesh, meshsize, "incircuit",
    # magdirection, group, turns)
    # From manual: Set the selected block labels to have the properties:
    # – Block property "blockname".
    # – automesh: 0 = mesher defers to mesh size constraint defined in
    # meshsize, 1 = mesher automatically chooses the mesh density.
    # – meshsize: size constraint on the mesh in the block marked by this label
    # – Block is a member of the circuit named "incircuit"
    # – The magnetization is directed along an angle in measured in degrees
    # denoted by the parameter magdirection. Alternatively, magdirection can be
    # a string containing a formula that prescribes the magnetization direction
    # as a function of element position. In this formula theta and R denotes
    # the angle in degrees of a line connecting the center each element with
    # the origin and the length of this line, respectively; x and y denote the
    # x- and y-position of the center of the each element. For axisymmetric
    # problems, r and z should be used in place of x and y.
    # – A member of group number group
    # – The number of turns associated with this label is denoted by turns.

    femm.mi_selectlabel(*etiquettes_dict['pcb_prim'])
    femm.mi_selectlabel(*etiquettes_dict['pcb_sec'])
    femm.mi_setblockprop('fr4', 1, 0, '<None>', 0, 1, 0)
    femm.mi_clearselected()

    femm.mi_selectlabel(*etiquettes_dict['air'])
    femm.mi_setblockprop('air', 1, 0, '<None>', 0, 1, 0)
    femm.mi_clearselected()

    femm.mi_selectlabel(*etiquettes_dict['gel'])
    femm.mi_setblockprop('silgel', 1, 0, 'None')
    femm.mi_clearselected()

    femm.mi_selectlabel(*etiquettes_dict['isolant'])
    femm.mi_setblockprop(tg.material_dielectric, 1, 0, 'None')
    femm.mi_clearselected()


def calc_inductance(tg, currents, inductances=(None, None), **kwargs):
    ''' Setup of magneto-static problem in femm to calculate inductance and
    resistance of planar transformer.

    Args:
        tg (:obj:'TransformerGeometry'): tg contains all geometry information
        of the transformer.
        currents (list of float): currents in the primary and secondary side
        circuits on the form [I_prim, I_sec].
        inductances (list of float): self-inductances of the primary and
        secondary side circuits on the form [L_prim, L_sec]. Used to calculate
        mutual inductance.

    Returns:
        (inductance, resistance): calculated self or mutual inductance and
        equivalent series resistance of either primary or secondary circuit.
    '''

    etiquettes_dict = {}  # Dictionary to store coordinates of nodes

    # initialitiation of the magneto-static problem
    boundary_radius = 2 * tg.radius_dielectric
    initial_setup(boundary_radius, currents, **kwargs)

    # draw geometry and add block labels
    add_conductors(tg, etiquettes_dict)
    add_pcbs(tg)
    add_isolation(tg)
    add_block_labels(tg, etiquettes_dict)

    # mi zoomnatural()
    # From manual: zooms to a “natural” view with sensible extents.
    femm.mi_zoomnatural()

    # Saving geometry file
    femm.mi_saveas('inductance_transformer.fem')

    # Meshing and analysis
    # From manual: Note that it is not necessary to run mesh before performing
    # an analysis, as mi_analyze() will make sure the mesh is up to date before
    # running an analysis.

    # mi analyze(flag)
    # From manual: runs fkern to solve the problem. The flag parameter controls
    # whether the fkern window is visible or minimized. For a visible window,
    # either specify no value for flag or specify 0. For a minimized window,
    # flag should be set to 1.
    femm.mi_analyze(1)

    # Post-processing
    femm.mi_loadsolution()

    # mo_seteditmode(mode)
    # From manual: Sets themode of the postprocessor to point, contour, or area
    # mode. Valid entries for mode are "point", "contour", and "area".
    femm.mo_seteditmode('area')

    # mo_blockintegral(type)
    # From manual: Calculate a block integral for the selected blocks
    # Type Definition
    # 0 A · J
    # 1 A
    # 2 Magnetic field energy
    # 3 Hysteresis and/or lamination losses
    # 4 Resistive losses
    # 5 Block cross-section area
    # 6 Total losses
    # 7 Total current
    # 8 Integral of Bx (or Br) over block
    # 9 Integral of By (or rBz) over block
    # 10 Block volume
    # ...

    # mo_getcircuitproperties("circuit")
    # From manual: Used primarily to obtain impedance information associated
    # with circuit properties. Properties are returned for the circuit property
    # named "circuit". Three values are returned by the function. In order,
    # these results are:
    # – current Current carried by the circuit
    # – volts Voltage drop across the circuit
    # – flux_re Circuit’s flux linkage

    # mo_groupselectblock(n)
    # From manual: Selects all the blocks that are labeled by block labels
    # that are members of group n. If no number is specified (i.e.
    # mo_groupselectblock() ), all blocks are selected.

    # Calculate the inductance of the circuit with non-zero current. If both
    # currents are given, we calculate the mutual inductance.
    L1, L2 = inductances
    if (currents[0] > 0) and (currents[1] == 0):
        circ = femm.mo_getcircuitproperties('phase_prim')
        resistance = circ[1].real
        inductance = abs(circ[2] / circ[0])
    elif (currents[0] == 0) and (currents[1] > 0):
        circ = femm.mo_getcircuitproperties('phase_sec')
        resistance = circ[1].real
        inductance = abs(circ[2] / circ[0])
    else:
        femm.mo_groupselectblock()
        # axisymmetric problem, integral is multiplied by 2
        Wm = femm.mo_blockintegral(2) * 2
        inductance = ((Wm - 0.5 * (L1*currents[1]**2 + L2*currents[0]**2)) /
                      (currents[0] * currents[1]))
        resistance = 0
        femm.mo_clearblock()

    if kwargs.get('close') is True:
        femm.closefemm()

    return (inductance, resistance)
