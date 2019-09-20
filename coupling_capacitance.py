# -*- coding: utf-8 -*-
"""
Created on Wed Dec 05 15:30:49 2018

@author: olechrs

###########################################################################
# Simulation FEMM lancé par Matlab. Il faut spécifier la géometrie
# des spires (16 valeur) du transformateur où le primaire et le secondaire
# sont fabriqués comme des circuits imprimés, geoemtrie[16].


        Matrices: index 1 pour le primaire (en haut),
                    index 2 pour le secondaire (en bas)

        |                    l_cuivre        l_cuivre
        |                   <----> l_entre <---->
        |<--r_interieur(1)-> ____ <-------> ____        no_sp_prim == 2
        |___________________|____|_________|____|_____
 y=0___|_____________________________________________| ep_pcb
        |                           |    |  /\
        |                            |____|  \/ep_cuivre
        |<-----r_interieur(2)------>                no_sp_sec == 1
        |
        |<---------------r_pcb------------------------>
        |
 l'axe de symétrie

###########################################################################
"""

import femm
import numpy as np


def initial_setup(limite_externe):
    # Commandes liés au logiciel pour acceder toutes les commandes FEMM
    # Les commandes prochaines sont lancé par le fichier principal de
    # l'algorithme génétique pour eviter des éxecuter plusieurs fois. En plus,
    # le fichier openfemm.m est modifié pour éviter que la fenêtre FEMM
    # travaille à l'arrière-plan: le drapeu '-windowhide' est ajouté à la
    # commande au système.
    femm.openfemm(1)
    # mi_minimize()
    # From manual: minimizes the active magnetics input view.
    femm.main_minimize()

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
    precision = 1e-9
    femm.ei_probdef('millimeters', 'axi', precision, 100, 30)

    # Circuit
    # From manual: ei_addconductorprop("conductorname", Vc, qc, conductortype)
    # adds a new conductor property with name "conductorname" with either a
    # prescribed voltage or a prescribed total charge. Set the unused property
    # to zero. The conductortype parameter is 0 for prescribed charge and 1 for
    # prescribed voltage
    femm.ei_addconductorprop('one', 1, 0, 1)
    femm.ei_addconductorprop('zero', 0, 0, 1)

    # Trace de la geometrie
    # ei_drawrectangle(x1, y1, x2, y2)
    # From manual: no discription

    # ei selectsegment(x,y)
    # From manual: Select the line segment closest to (x,y)

    # ei_setsegmentprop("propname", elementsize, automesh, hide,
    #                   group, "inconductor",)
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
    # ei addmaterial("materialname", er er ?)
    # TODO fill in missing description
    femm.ei_addmaterial('air', 1, 1, 0)
    femm.ei_addmaterial('FR4', 4.4, 4.4, 0)
    femm.ei_addmaterial('Polysterimide', 3.5, 3.5, 0)
    femm.ei_addmaterial('Teflon', 2.1, 2.1, 0)
    femm.ei_addmaterial('Silgel', 2.7, 2.7, 0)

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
    femm.ei_makeABC(7, limite_externe, 0, 0, 0)


def translate_geometry(g):
    ''' Translates from old to new style of handling geometry parameters '''
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


def set_conductor_boundry(coords, volt):
    '''
    Set voltage of the four sides of a square conductor to either 0 or 1 V
    INPUT: corner, length and height of conductor
    '''
    # ei_selectsegment(x,y)
    # From manual: Select the line segment closest to (x,y)
    femm.ei_selectsegment(coords[0], np.average((coords[1], coords[3])))
    femm.ei_selectsegment(np.average((coords[0], coords[2])), coords[1])
    femm.ei_selectsegment(coords[2], np.average((coords[1], coords[3])))
    femm.ei_selectsegment(np.average((coords[0], coords[2])), coords[3])
    if volt == 1:
        femm.ei_setsegmentprop('<None>', 0, 1, 0, 0, 'one')
    elif volt == 0:
        femm.ei_setsegmentprop('<None>', 0, 1, 0, 0, 'zero')
    femm.ei_clearselected()


def set_conductor_label(coords, label_name, label_dict):
    '''
    Set a block label inside a square conductor, and give the block label
    properties. INPUT: corner, length and height of conductor, label name and
    label dictionary
    '''
    # ei_addblocklabel(x,y)
    # From manual: Add a new block label at (x,y)
    label_coord = (np.average((coords[0], coords[2])),
                   np.average((coords[1], coords[3])))
    femm.ei_addblocklabel(*label_coord)
    label_dict[label_name] = label_coord
    return label_dict


def draw_conductor(coords, in_conductor, label_name, label_dict):
    ''' Draw rectabngualr shapeof conducture in problem and add label '''
    femm.ei_drawrectangle(*coords)
    # ajouter la tension pour les trace du bobinage du primaire
    set_conductor_boundry(coords, in_conductor)
    # label du bobinage du primaire
    label_dict = set_conductor_label(coords, label_name, label_dict)
    return label_dict


def add_conductors(geometry, z0, dz, label_dict):
    ''' Add all conductors of a pcb to the problem '''
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
        height = np.cumsum([initial_value] + list(d_height)).tolist()
    else:
        print ('get_height_between_layers(): invalid number of layers')
    return height


def coords_of_rectangle(x0, y0, dx, dy):
    '''
    Return array with the coordinates of the two opposite corners of a
    rectangle. Useful as femm commands don't accept arrays as input.
    '''
    return (x0, y0, x0+dx, y0+dy)


def calc_coupling_capacitance(transformer_geometry):

    # Initiation de parametres
    geometry = translate_geometry(transformer_geometry)
    # Grandeurs de la géometrie
    no_spires_prim = geometry[0]
    no_couches_prim = geometry[1]
    no_spires_sec = geometry[2]
    no_couches_sec = geometry[3]
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

    initial_setup(limits_externes)

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

    femm.ei_drawrectangle(*coords_disc_dielectrique)
    femm.ei_drawrectangle(*coords_gel_dielectrique)
    femm.ei_drawrectangle(*coords_pcb_prim)
    femm.ei_drawrectangle(*coords_pcb_sec)

    # mi zoomnatural()
    # From manual: zooms to a “natural” view with sensible extents.
    femm.ei_zoomnatural

    # Ajoute de block labels, étiquettes dans les surfaces
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

    # label pour le PCB
    coords = [2, ep_dielectric + ep_pcb_noyau / 2.]
    femm.ei_addblocklabel(*coords)
    etiquettes_dict['pcb_prim'] = coords

    coords = [2, - ep_pcb_noyau / 2.]
    femm.ei_addblocklabel(*coords)
    etiquettes_dict['pcb_sec'] = coords

    # label pour l'air autour
    coords = [2, ep_dielectric * 2 + ep_gel]
    femm.ei_addblocklabel(*coords)
    etiquettes_dict['air'] = coords

    # label pour le gel autour
    coords = [r_pcb + 2, ep_dielectric + ep_gel / 2.]
    femm.ei_addblocklabel(*coords)
    etiquettes_dict['gel'] = coords

    # label pour l'isolant
    coords = [2, ep_dielectric / 2.]
    femm.ei_addblocklabel(*coords)
    etiquettes_dict['isolant'] = coords

    # Associer blocks avec materiaux
    # ei setblockprop("blockname", automesh, meshsize, group) Set the selected
    # block labels to have the properties: Block property "blockname".
    # automesh: 0 = mesher defers to mesh size constraint defined in meshsize,
    # 1 = mesher automatically chooses the mesh density. meshsize: size
    # constraint on the mesh in the block marked by this label. A member of
    # group number group

    femm.ei_selectlabel(*etiquettes_dict['pcb_prim'])
    femm.ei_selectlabel(*etiquettes_dict['pcb_sec'])
    femm.ei_setblockprop('FR4', 1, 0, 'None')
    femm.ei_clearselected()

    femm.ei_selectlabel(*etiquettes_dict['air'])
    femm.ei_setblockprop('air', 1, 0, 'None')
    femm.ei_clearselected()

    femm.ei_selectlabel(*etiquettes_dict['gel'])
    femm.ei_setblockprop('Silgel', 1, 0, 'None')
    femm.ei_clearselected()

    femm.ei_selectlabel(*etiquettes_dict['isolant'])
    if materiel_dielectrique == 1:
        femm.ei_setblockprop('FR4', 1, 0, 'None')
    elif materiel_dielectrique == 2:
        femm.ei_setblockprop('Polysterimide', 1, 0, 'None')
    elif materiel_dielectrique == 3:
        femm.ei_setblockprop('Teflon', 1, 0, 'None')
    femm.ei_clearselected()

    # materiau du bobinage du primaire
    for j in range(no_couches_prim):
        for i in range(no_spires_prim):
            identifier = no_spires_prim * j + i
            femm.ei_selectlabel(*etiquettes_dict['prim' + str(identifier)])
            femm.ei_setblockprop('<No Mesh>', 1, 0, 'phase_prim')
            femm.ei_clearselected()

    # materiau du bobinage du secondaire
    for j in range(no_couches_sec):
        for i in range(no_spires_sec):
            identifier = no_spires_sec * j + i
            femm.ei_selectlabel(*etiquettes_dict['sec' + str(identifier)])
            femm.ei_setblockprop('<No Mesh>', 1, 0, 'phase_sec')
            femm.ei_clearselected()

    # Enregistrement
    femm.ei_saveas('capacitance_transfo.fee')

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
    femm.ei_analyze(1)

    # Post-processeur
    femm.ei_loadsolution()

    # eo_seteditmode(mode)
    # From manual: Sets themode of the postprocessor to point, contour, or area
    # mode. Valid entries for mode are "point", "contour", and "area".
    femm.eo_seteditmode('area')

    # eo_selectblock(x,y)
    # From manual: Select the block that contains point (x,y).
    femm.eo_selectblock(*etiquettes_dict['pcb_prim'])
    femm.eo_selectblock(*etiquettes_dict['pcb_sec'])
    femm.eo_selectblock(*etiquettes_dict['air'])
    femm.eo_selectblock(*etiquettes_dict['gel'])
    femm.eo_selectblock(*etiquettes_dict['isolant'])

    # eo_blockintegral(type)
    # From manual: Calculate a block integral for the selected blocks
    # type Integral
    # 0 Stored Energy
    # 1 Block Cross-section
    # 2 Block Volume
    # 3 Average D over the block
    # 4 Average E over the block
    # 5 Weighted Stress Tensor Force
    # 6 Weighted Stress Tensor Torque
    # Returns one or two floating point values as results, e.g.:
    # Fx, Fy = eo blockintegral(4)

    # eo_getconductorproperties("conductor")Properties are returned for the
    # conductor property named ”conductor”. Two values are returned: The
    # voltage of the specified conductor, and the charge carried on the
    # specified conductor.

    # Calculer la capacité
#    energy = femm.eo_blockintegral(0)
    femm.eo_clearblock()
    circuit_properties = femm.eo_getconductorproperties('one')

#    capacitance = (2 * energy[0], circuit_properties[1])
    capacitance = circuit_properties[1]

    # Commandes liés au logiciel pour fermer correctement le FEMM
    # femm.eo_close()
    # femm.ei_close()
    return capacitance
