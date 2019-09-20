# -*- coding: utf-8 -*-
"""
Created on Wed Dec 05 15:30:49 2018

@author: olechrs

###########################################################################
# Simulation FEMM lancé par Matlab. Il faut spécifier la géometrie
# des spires (16 valeur) du transformateur où le primaire et le secondaire
# sont fabriqués comme des circuits imprimés, geoemtrie[16].
#
#
#        Matrices: index 1 pour le primaire (en haut),
#                    index 2 pour le secondaire (en bas)
#
#        |                    l_cuivre        l_cuivre
#        |                   <----> l_entre <---->
#        |<--r_interieur(1)-> ____ <-------> ____        no_sp_prim == 2
#        |___________________|____|_________|____|_____
# y=0___|_____________________________________________| ep_pcb
#        |                           |    |  /\
#        |                            |____|  \/ep_cuivre
#        |<-----r_interieur(2)------>                no_sp_sec == 1
#        |
#        |<---------------r_pcb------------------------>
#        |
# l'axe de symétrie

 Geometry array
 1 no_spires_prim =      nombre de spires au primaire
 2 no_couches_prim =     nombre de couche au primaire
 3 no_spires_sec =       nombre de spires au secondaire
 4 no_couches_sec =      nombre de couche au secondaire
 5 ep_pcb_noyau =        epaisseur noyau du circuit imprimé cuivre non-inclu
                         [mm]
 6 ep_pcb_prepreg =      epaisseur couches isolantes additionelles du circuit
                         imprimé [mm]
 7 ep_dielectric =       epaisseur de la couche isolante principale [mm]
 8 ep_cuivre =           epaisseur cuivre [mm]
 9 l_cuivre =            largeur de fil en cuivre [mm]
10 r_interieur =         radii entre l'axe et le premier tour [mm]
11 l_entre =             largeur entre chaque tour [mm]
12 r_pcb =               radii de pcb [mm]
13 r_dielectrique =      radii du disc isolant [mm]
14 materiel_dielectrique =  materiel_dielectrique entre le primaire et le
                         secondaire du transformateur. 'FR4', 'polyesterimide'
                         et 'teflon' sont acceptés
15 ep_gap =              gap entre le cuivre sur le PCB et l'isolant [mm]
16 ep_gel =              epaisseur du gel à chaque coté de l'isolant [mm]
17 r_gel =              radii du gel supérieure au disc isolant [mm]
###########################################################################
"""

import femm
import numpy as np


def initial_setup(limite_externe, courants):
    # Commandes liés au logiciel pour acceder toutes les commandes FEMM
    # Les commandes prochaines sont lancé par le fichier principal de
    # l'algorithme génétique pour eviter des éxecuter plusieurs fois. En plus,
    # le fichier openfemm.m est modifié pour éviter que la fenêtre FEMM
    # travaille à l'arrière-plan: le drapeu '-windowhide' est ajouté à la
    # commande au système.
    femm.openfemm()
    # mi_minimize()
    # From manual: minimizes the active magnetics input view.
    #femm.main_minimize()

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
    frequence = 6.78e6                # frequence de coupure
    precision = 5e-9                    # precision du calcul de la simulation
    femm.mi_probdef(frequence, 'millimeters', 'axi', precision, 50, 30)

    # Circuit
    # mi addcircprop("circuitname", i, circuittype)
    # From manual: adds a new circuit property with name "circuitname" with a
    # prescribed current, i. The circuittype parameter is 0 for a
    # parallel-connected circuit and 1 for a series-connected circuit.

    # Les courants donnés sont pour le primaire et le secondaire
    I1, I2 = courants

    femm.mi_addcircprop('phase_prim', I1, 1)
    femm.mi_addcircprop('phase_sec', I2, 1)

    # Trace de la geometrie
    # mi_drawrectangle(x1, y1, x2, y2)
    # From manual: no discription

    # ei selectsegment(x,y)
    # From manual: Select the line segment closest to (x,y)

    # mi_setsegmentprop("propname",
    #                   elementsize, automesh, hide, group, "inconductor",)
    # From manual: Set the select segments to have:
    # Boundary property "propname"
    # Local element size along segment no greater than elementsize
    # automesh: 0 = mesher defers to the element constraint defined by
    # elementsize, 1 = mesher automatically chooses mesh size along the
    # selected segments
    # hide: 0 = not hidden in post-processor, 1 == hidden in post processor
    # A member of group number group
    # A member of the conductor specified by the string "inconductor". If the
    # segment is not part of a conductor, this parameter can be specified as
    # "<None>".

    # Materiaux
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
    femm.mi_addmaterial('FR4', 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    femm.mi_addmaterial('cuivre', 1, 1, 0, 0, 58, 0, 0, 0, 0, 0, 0)
    femm.mi_addmaterial('Polysterimide', 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    femm.mi_addmaterial('Teflon', 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    femm.mi_addmaterial('Silgel', 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)

    # Conditions limites
    # ei makeABC(n,R,x,y,bc)
    # From manual: creates a series of circular shells that emulate the
    # impedance of an unbounded domain (i.e. an Improvised Asymptotic Boundary
    # Condition). The n parameter contains the number of shells to be used
    # (should be between 1 and 10), R is the radius of the solution domain, and
    # (x,y) denotes the center of the solution domain. The bc parameter should
    # be specified as 0 for a Dirichlet outer edge or 1 for a Neumann outer
    # edge. If the function is called without all the parameters, the function
    # makes up reasonable values for the missing parameters.
    femm.mi_makeABC(7, limite_externe, 0, 0, 0)


def translate_geometry(g):
    ''' Translates from old to new style of handling geometry parameters'''
    geometry = []
    geometry.append(g.turns_primary)
    geometry.append(g.layers_primary)
    geometry.append(g.turns_secondary)
    geometry.append(g.layers_secondary)
    geometry.append(g.height_pcb_core)
    geometry.append(g.height_pcb_prepreg)
    geometry.append(g.height_dielectric)
    geometry.append(g.height_copper)
    geometry.append(g.width_copper)
    geometry.append(g.radius_inner_track)
    geometry.append(g.width_between_tracks)
    geometry.append(g.radius_pcb)
    geometry.append(g.radius_dielectric)
    if g.material_dielectric == 'fr4':
        geometry.append(1)
    elif g.material_dielectric == 'polysterimide':
        geometry.append(2)
    elif g.material_dielectric == 'teflon':
        geometry.append(3)
    else:
        geometry.append(-1)
    geometry.append(g.height_gap)
    geometry.append(g.height_gel)
    geometry.append(g.radius_gel)

    return geometry


def coords_of_rectangle(x0, y0, dx, dy):
    '''
    Return array with the coordinates of the two opposite corners of a
    rectangle. Useful as femm commands don't accept arrays as input.
    '''
    return (x0, y0, x0+dx, y0+dy)


def get_height_between_layers(no_layers, initial_value, d_height):
    '''
    return an array with the height between copper layers in the PCB depending
    on the number of layers and geometry parameters
    '''
    if no_layers == 1:
        height = (initial_value,)
    elif no_layers == 2:
        height = (initial_value, initial_value + d_height)
    elif no_layers == 4:
        height = np.cumsum([initial_value]+list(d_height)).tolist()
    else:
        print ('get_height_between_layers(): invalid number of layers')
    return height


def set_conductor_label(coords, in_conductor, label_name, label_dict):
    '''
    set a block label inside a square conductor, and give the block label
    properties. INPUT: corner, length and height of conductor, label name and
    label dictionary
    '''
    # mi_addblocklabel(x,y)
    # From manual: Add a new block label at (x,y)
    label_coord = (np.average((coords[0], coords[2])), np.average((coords[1],
                                                                   coords[3])))
    femm.mi_addblocklabel(*label_coord)
    femm.mi_selectlabel(*label_coord)
    if in_conductor == 1:
        femm.mi_setblockprop('cuivre', 1,  0, 'phase_prim', 0, 1, 1)
    elif in_conductor == 0:
        femm.mi_setblockprop('cuivre', 1,  0, 'phase_sec', 0, 1, 1)
    femm.mi_clearselected()
    label_dict[label_name] = label_coord
    return label_dict


def draw_conductor(coords, in_conductor, label_name, label_dict):
    ''' '''
    femm.mi_drawrectangle(*coords)
    # ajouter la tension pour les trace du bobinage du primaire
    # label du bobinage du primaire
    label_dict = set_conductor_label(coords,
                                     in_conductor,
                                     label_name,
                                     label_dict)
    return label_dict


def add_conductors(geometry, z0, dz, label_dict):
    '''  Add conductors on a pcb layer to the problem '''
    no_spires_prim = geometry[0]
    no_couches_prim = geometry[1]
    no_spires_sec = geometry[2]
    no_couches_sec = geometry[3]
    ep_cuivre = geometry[7]
    l_cuivre = geometry[8]
    r_interieur = geometry[9]
    l_entre = geometry[10]

    ep_couches = get_height_between_layers(no_couches_prim, z0[0], dz)
    for c in range(no_couches_prim):
        for i in range(no_spires_prim):
            # coordiantes of conductor
            coords_conductor = coords_of_rectangle(r_interieur[0] +
                                                   (l_entre[0] + l_cuivre) * i,
                                                   ep_couches[c], l_cuivre,
                                                   ep_cuivre)
            # desiner la conductrice, ajouter la tension  et le label
            label_dict = draw_conductor(coords_conductor, 1,
                                        'prim' + str(i + c * no_spires_prim),
                                        label_dict)
    # définition du bobinage du secondaire
    ep_couches = get_height_between_layers(no_couches_prim, z0[1], -dz)
    for c in range(no_couches_sec):
        for i in range(no_spires_sec):
            coords_conductor = coords_of_rectangle(r_interieur[1] +
                                                   (l_entre[1] + l_cuivre) * i,
                                                   ep_couches[c], l_cuivre,
                                                   -ep_cuivre)
            # desiner la conductrice, ajouter la tension  et le label
            label_dict = draw_conductor(coords_conductor, 0,
                                        'sec' + str(i + c * no_spires_sec),
                                        label_dict)


def calc_inductance(transformer_geometry, courants, inductances=(None, None)):

    # Initiation de parametres
    # Si les courants sont plus grand que zero et les inductances du primaire
    # et du secondaire sont données, on calcul l'inductance mutuel.
    L1, L2 = inductances
    L, R = [0, 0]
    geometry = translate_geometry(transformer_geometry)

    # Grandeurs de la géometrie
#    no_spires_prim = geometry[0]
    no_couches_prim = geometry[1]
#    no_spires_sec = geometry[2]
#    no_couches_sec = geometry[3]
    ep_pcb_noyau = geometry[4]
    ep_pcb_prepreg = geometry[5]
    ep_dielectric = geometry[6]
    ep_cuivre = geometry[7]
#    l_cuivre = geometry[8]
#    r_interieur = geometry[9]
#    l_entre = geometry[10]
    r_pcb = geometry[11]
    r_dielectrique = geometry[12]
    materiel_dielectrique = geometry[13]
    ep_gap = geometry[14]
    ep_gel = geometry[15]
    r_gel = geometry[16]

    # condition aux limits, limit externe de la simulation
    limits_externes = 2 * r_dielectrique

    # dictionnaire pour gérer les etiquette de la géometrie
    etiquettes_dict = {}

    # initialitiation de la problème
    initial_setup(limits_externes, courants)

    # définition de l'dielectrique
    coords_disc_dielectrique = coords_of_rectangle(0, 0, r_dielectrique,
                                                   ep_dielectric)
    # définition du gel dielectrique
    coords_gel_dielectrique = coords_of_rectangle(0, -ep_gel, r_gel,
                                                  2*ep_gel+ep_dielectric)
    # définition du PCB
    coords_pcb_prim = 0
    coords_pcb_sec = 0
    z0 = (ep_dielectric + ep_gap, -ep_gap)

    # dessiner le bobinage
    if no_couches_prim == 1:
        # print('une couche par PCB')
        # définition du bobinage
        dz = 0
        add_conductors(geometry, z0, dz, etiquettes_dict)
        # définition du PCB au primaire
        coords_pcb_prim = coords_of_rectangle(0, ep_dielectric + ep_cuivre +
                                              ep_gap, r_pcb, ep_pcb_noyau)
        # définition du PCB au secondaire
        coords_pcb_sec = coords_of_rectangle(0, -ep_cuivre - ep_gap, r_pcb,
                                             -ep_pcb_noyau)

    elif no_couches_prim == 2:
        # print('deux couches par PCB')
        # définition du bobinage
        dz = ep_cuivre + ep_pcb_noyau
        add_conductors(geometry, z0, dz, etiquettes_dict)
        # définition du PCB au primaire
        coords_pcb_prim = coords_of_rectangle(0, ep_dielectric + ep_cuivre +
                                              ep_gap, r_pcb, ep_pcb_noyau)
        # définition du PCB au secondaire
        coords_pcb_sec = coords_of_rectangle(0, -ep_cuivre - ep_gap, r_pcb,
                                             -ep_pcb_noyau)

    elif no_couches_prim == 4:
        # print('quatre couches par PCB')
        # définition du bobinage
        dz = np.array((ep_cuivre + ep_pcb_prepreg, ep_cuivre + ep_pcb_noyau,
                       ep_cuivre + ep_pcb_prepreg))
        add_conductors(geometry, z0, dz, etiquettes_dict)
        # définition du isolant au primaire
        coords_pcb_prim = coords_of_rectangle(0, ep_dielectric + ep_cuivre +
                                              ep_gap, r_pcb, ep_pcb_noyau +
                                              2*ep_pcb_prepreg + 3*ep_cuivre)
        # définition du isolant au secondaire
        coords_pcb_sec = coords_of_rectangle(0, -ep_cuivre - ep_gap, r_pcb,
                                             -ep_pcb_noyau - 2*ep_pcb_prepreg -
                                             3*ep_cuivre)
    else:
        'paramètre de couche inconnu'
        return 0

    femm.mi_drawrectangle(*coords_disc_dielectrique)
    femm.mi_drawrectangle(*coords_gel_dielectrique)
    femm.mi_drawrectangle(*coords_pcb_prim)
    femm.mi_drawrectangle(*coords_pcb_sec)

    # mi zoomnatural()
    # From manual: zooms to a “natural” view with sensible extents.
    femm.mi_zoomnatural

    # Ajoute de block labels, étiquettes dans les surfaces
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

    # label pour le PCB
    coords = [2, ep_dielectric + ep_pcb_noyau / 2.]
    femm.mi_addblocklabel(*coords)
    etiquettes_dict['pcb_prim'] = coords

    coords = [2, - ep_pcb_noyau / 2.]
    femm.mi_addblocklabel(*coords)
    etiquettes_dict['pcb_sec'] = coords

    # label pour l'air autour
    coords = [2, ep_dielectric * 2 + ep_gel]
    femm.mi_addblocklabel(*coords)
    etiquettes_dict['air'] = coords

    # label pour le gel autour
    coords = [r_pcb + 2, ep_dielectric + ep_gel / 2.]
    femm.mi_addblocklabel(*coords)
    etiquettes_dict['gel'] = coords

    # label pour l'isolant
    coords = [2, ep_dielectric / 2.]
    femm.mi_addblocklabel(*coords)
    etiquettes_dict['isolant'] = coords

    # Associer blocks avec materiaux
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
    femm.mi_setblockprop('FR4', 1, 0, '<None>', 0, 1, 0)
    femm.mi_clearselected()

    femm.mi_selectlabel(*etiquettes_dict['air'])
    femm.mi_setblockprop('air', 1, 0, '<None>', 0, 1, 0)
    femm.mi_clearselected()

    femm.mi_selectlabel(*etiquettes_dict['gel'])
    femm.mi_setblockprop('Silgel', 1, 0, 'None')
    femm.mi_clearselected()

    femm.mi_selectlabel(*etiquettes_dict['isolant'])
    if materiel_dielectrique == 1:
        femm.mi_setblockprop('FR4', 1, 0, 'None')
    elif materiel_dielectrique == 2:
        femm.mi_setblockprop('Polysterimide', 1, 0, 'None')
    elif materiel_dielectrique == 3:
        femm.mi_setblockprop('Teflon', 1, 0, 'None')
    femm.mi_clearselected()

    # Enregistrement
    femm.mi_saveas('inductance_transfo.fem')

    # Maillage
    # From manual: Note that it is not necessary to run mesh before performing
    # an analysis, as mi_analyze() will make sure the mesh is up to date before
    # running an analysis.

    # Resolution
    # mi analyze(flag)
    # From manual: runs fkern to solve the problem. The flag parameter controls
    # whether the fkern window is visible or minimized. For a visible window,
    # either specify no value for flag or specify 0. For a minimized window,
    # flag should be set to 1.
    femm.mi_analyze(1)

    # Post-processeur
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

    # Calcule les parametres du primaire
    if (courants[0] > 0) and (courants[1] == 0):
        # Power_losses = abs(femm.mo_blockintegral(6))
        # R = Power_losses / I1**2
        circ = femm.mo_getcircuitproperties('phase_prim')
        R = circ[1].real
        L = abs(circ[2] / circ[0])
    # Calcule les parametres du secondaire
    elif (courants[0] == 0) and (courants[1] > 0):
        # Power_losses = abs(femm.mo_blockintegral(6))
        # R = Power_losses / I2**2
        circ = femm.mo_getcircuitproperties('phase_sec')
        R = circ[1].real
        L = abs(circ[2] / circ[0])
#        print (circ)
    # Calcule les parametres mutuelle
    else:
        femm.mo_groupselectblock()
        # axisymmetric problem, integral is multiplied by 2
        Wm = femm.mo_blockintegral(2) * 2
        L = (Wm - 0.5*(L1*courants[1]**2 + L2*courants[0]**2)) / (courants[0] *
                                                                  courants[1])
        R = 0
    femm.mo_clearblock()

    # mo_showdensityplot(legend,gscale,upper_B,lower_B,type)
    # From manual: Shows the flux density plot with options:
    # – legend Set to 0 to hide the plot legend or 1 to show the plot legend.
    # – gscale Set to 0 for a colour density plot or 1 for a grey scale 
    #   density plot.
    # – upper_B Sets the upper display limit for the density plot.
    # – lower_B Sets the lower display limit for the density plot.
    # – type Type of density plot to display. Valid entries are ’mag’, ’real’,
    # and ’imag’ for magnitude, real component, and imaginary component of B,
    # respectively. Alternatively, current density can be displayed by
    # specifying ’jmag’, ’jreal’, and ’jimag’ for magnitude, real component,
    # and imaginary component of J, respectively.

#    femm.mo_showdensityplot(1, 0, 2e-3, 0, 'mag')
    # Commandes liés au logiciel pour fermer correctement le FEMM
#    femm.mo_close()
#    femm.mi_close()
#    femm.closefemm()

    return (L, R)
