#!/usr/bin/env python3

from vaspwfc import vaspwfc
from vasp_constant import *
import numpy as np
import subprocess
from aflow_sym_python import Symmetry
import json
import bz2
import math as m

def get_pointgroup(pos_file, settings):
    """
    Gets the point group of a crystal structure via AFLOW-SYM.

    Input:
        pos_file: string that is the path to a crystal structure file like POSCAR or CONTCAR
        settings: loaded json settings file

    Returns:
        Point group name as string
    """


    #pg = subprocess.check_output("aflow --aflowsym="+str(settings['aflow_tolerance'])+" < "+pos_file+" | awk '/Schoenflies/ {print $NF}' | tail -1", shell=True, text=True)
    pg = subprocess.run("aflow --pgroup_xtal="+str(settings['aflow_tolerance'])+ \
    " < "+pos_file+" | awk '/Schoenflies/ {print $NF}' | tail -1", shell=True, \
    check=True, text=True, capture_output=True).stdout
    subprocess.run("rm aflow.*.out.xz" ,shell=True)

    pg = group_name_conv(pg)

    return pg

def get_symmetry_operators_old(sym, pos_file, settings):
    """
    Fetches symmetry operator names/symbols and matrices

    Input:
        sym: An object of the Symmetry class (from aflow_sym_python)
        pos_file: string that is the path to a crystal structure file like POSCAR or CONTCAR
        settings:

    Returns:
        Schoenflies symbol of transformation
        Transformation matrix
        Axis of rotation
        Angle rotated
    """

    out = sym.get_symmetry(pos_file,settings['aflow_tolerance'])['pgroup_xtal']



    l = []
    for i in range(len(out)):
        l.append((out[i]['Schoenflies'],i))
    l.sort(key=lambda tup: tup[0])

    ran = []
    for sym, i in l:
        ran.append(i)

    symbols = []
    matrices = []
    angle = []
    axis = []

    for i in ran:#range(len(out)):
        symbols.append(out[i]['Schoenflies'])
        matrices.append(out[i]['Uf'])
        angle.append(out[i]['angle'])
        axis.append(out[i]['axis'])
    return symbols, matrices, axis, angle

def get_symmetry_operators(sym, pos_file, PGname, settings):
    """
    Fetches symmetry operator names/symbols and matrices, and
    sorts them into conjugacy classes like in the character table.

    Input:
        sym: An object of the Symmetry class (from aflow_sym_python)
        pos_file: string that is the path to a crystal structure file like POSCAR or CONTCAR
        PGname: name of point group
        settings:

    Returns:
        Schoenflies symbol of transformation
        Transformation matrix
        Axis of rotation
        Angle rotated
    """

    out_json = subprocess.check_output("aflow --pgroup_xtal="+str(settings['aflow_tolerance'])+" --print=json --screen_only < "+pos_file, shell=True)

    out = json.loads(out_json)['pgroup_xtal']

    l = []
    for i in range(len(out)):
        l.append((out[i]['Schoenflies'],i))
    l.sort(key=lambda tup: tup[0])

    ran = []
    for sym, i in l:
        ran.append(i)

    symbols1 = []
    matrices1 = []
    angle1 = []
    axis1 = []

    for i in ran:
        symbols1.append(out[i]['Schoenflies'])
        matrices1.append(out[i]['Uf'])
        angle1.append(out[i]['angle'])
        axis1.append(out[i]['axis'])

    ch_table, pos_vec_rep = get_character_table(PGname,settings)
    class_symbols = ch_table[0][1:]
    mult_right_order = []
    classes_right_order = []
    int_previous = False
    for foo in class_symbols:
        try:
            mult_right_order.append(int(foo))
            int_previous = True
        except:
            classes_right_order.append(foo)
            if not int_previous:
                mult_right_order.append(1)
            int_previous = False

    class_symbols, class_perm, mult = sort_into_classes(matrices1, symbols1)
    perm = []

    for i, symb in enumerate(classes_right_order):
        symb_2 = "X"
        if symb == "E":
            symb = "1"
        elif symb == "sv" or symb == "sd" or symb == "sh":
            symb = "S2"
            symb_2 = "s"
        elif symb == "S6": # !!!!!????? S_n vs S_2n
            symb = "S3"
        elif symb == "C2'" or symb == "C2''" or symb == 'C2"':
            symb = "C2"

        occurance = [j for j, n in enumerate(class_symbols) if n == symb or n == symb_2]
        for ind in occurance:
            if mult_right_order[i] == mult[ind] and not ind in perm:
                perm.append(ind)

    class_perm2 = []
    for p in perm:
        class_perm2.append(class_perm[p])
    permutation = []
    for cp2 in class_perm2:
        for p in cp2:
            permutation.append(p)
    symbols = []
    matrices = []
    angle = []
    axis = []

    for i in permutation:
        symbols.append(symbols1[i])
        matrices.append(matrices1[i])
        angle.append(angle1[i])
        axis.append(axis1[i])

    return symbols, matrices, axis, angle

def sort_into_classes(sym_ops, symbols):
    """
    Sorts symmetry operators into conjugacy classes,

    Input:
        sym_ops: array of symmetry operators (3x3 matrix)
        symbols: array of symbols of each operator (e.g. C3, S2, i)
    Returns:
        class_symbols: array with symbols of each class
        class_perm: array with permutation of the operators into the classes
        mult: array with multiplicity (order) of each class
    """
    class_symbols = []
    class_perm = []
    used = []

    for i, A in enumerate(sym_ops):
        if i not in used:
            cls_perm = [i]
            class_symbols.append(symbols[i])
            used.append(i)
            for j, B in enumerate(sym_ops):
                for X in sym_ops:
                    conj = np.linalg.inv(X).dot(np.array(A).dot(X))
                    if np.all(B == conj) and j not in used:
                        cls_perm.append(j)
                        used.append(j)
            class_perm.append(cls_perm)

    mult = []
    for i in class_perm:
        mult.append(len(i))

    return class_symbols, class_perm, mult

def get_fixed_atoms(pos_file):
    """
    Calculates the location of the atoms that are fixed (mapped to themselves,
    possibly in an adjacent unit cell) for all symmetry transformations of a point group.

    Input:
        pos_file: string that is the path to a crystal structure file like POSCAR or CONTCAR

    Returns:
        Array with coordinates of fixed atoms
    """
    perms = []
    sym = Symmetry()
    fgroup = sym.get_symmetry(pos_file)['fgroup']

    fixed_point = [0,0,0]

    map1 = fgroup[1]['basis_atoms_map']
    natoms = len(map1)
    fixed_atoms = []

    if len(fgroup) > 0:
        for i, n in enumerate(map1):
            if n == i:
                fixed_atoms.append(i)

    print(fixed_atoms)

    to_be_removed = []

    if len(fgroup) > 2:
         for data in fgroup[2:]:
            map = data['basis_atoms_map']
            for n in fixed_atoms:
                if n != map[n] and n not in to_be_removed:
                    to_be_removed.append(n)

    for n in to_be_removed:
        fixed_atoms.remove(n)

    if len(fixed_atoms) == 0:
        print("No atoms at fixed points!")
        # Make function to check cycles of the permutation

    f = open(pos_file,"r")
    lines = f.readlines()[8:]

    try:
        float(lines[0].split()[0])
    except:
        lines = lines[1:]

    lines = lines[0:natoms-1]

    fixedpoints = []
    for i, n in enumerate(fixed_atoms):
        #print(lines[n])
        fp = lines[fixed_atoms[i]].split()[0:3]
        fixed_point = [float(fp[0]), float(fp[1]), float(fp[2])]
        fixedpoints.append(fixed_point)

    #fp = lines[fixed_atoms[0]].split()[0:3]
    #fixed_point = [float(fp[0]), float(fp[1]), float(fp[2])]

    return fixedpoints

def get_single_species(pos_file):
    """
    gets the coordinates of an atom of a unique species

    Input:
        pos_file: string that is the path to a crystal structure file like POSCAR or CONTCAR

    Returns:
        coordinates
    """
    f = open(pos_file,"r")
    lines = f.readlines()

    try:
        amounts = np.array(lines[6].split()).astype(np.float)
    except:
        return None

    if np.any(amounts == 1):
        ind = list(amounts).index(1)
        index = int(sum(amounts[0:ind]))
        temp = lines[8+index].split()[0:3]
        coord = [float(temp[0]), float(temp[1]), float(temp[2])]
        return coord

    return None

def find_average_position_general(realwf_array, procent):
    """
    Finds the average position of the wavefuntion (weighted center of mass)
    Input:
        realwf_array: array with realspace wavefuntions as output by vaspwfc.wfc_r
        procent: sets psi=0 when |psi|^2 < max * procent
    Returns:
        average position
    """
    sh = realwf_array[0].shape
    wf2_tot = np.zeros(sh)

    for wf in realwf_array:
        wf2_tot += wf*wf.conj()

    cutoff = np.max(wf2_tot)*procent

    [a1,a2,a3] = np.meshgrid(np.arange(sh[0])/(sh[0]-1),np.arange(sh[1])/(sh[1]-1),np.arange(sh[2])/(sh[2]-1), indexing="ij")

    wf2_cut = np.where(wf2_tot > cutoff, wf2_tot, 0)
    wf2_cut = wf2_cut/np.sum(wf2_cut)
    avg_pos = [np.sum(a1*wf2_cut), np.sum(a2*wf2_cut), np.sum(a3*wf2_cut)]

    return avg_pos

def find_circular_mean_realspace_opt(realwf_array, procent):
    """
    Finds the circular mean of the realspace wavefuntion squared, optially of low values of psi are cut off.
    Input:
        realwf: realspace wavefuntion as output by vaspwfc.wfc_r
        procent: sets psi=0 when |psi|^2 < max * procent
    Returns:
        position of the circular mean in the lattice basis
    """
    sh = realwf_array[0].shape

    wf2_tot = np.zeros(sh)

    for wf in realwf_array:
        wf2_tot += wf*wf.conj()

    cutoff = np.max(wf2_tot)*procent
    wf2_cut = np.where(wf2_tot > cutoff, wf2_tot, 0)
    wf2_cut = wf2_cut/np.sum(wf2_cut)

    grid = np.meshgrid(np.arange(sh[0])/(sh[0]-1),np.arange(sh[1])/(sh[1]-1),np.arange(sh[2])/(sh[2]-1), indexing="ij")

    circ_mean = []
    for a in grid:
        temp = np.angle(np.sum(m.e**(2j*m.pi*a)*wf2_cut))/(2*m.pi)
        if temp < 0:
            temp += 1
        circ_mean.append(temp)

    return circ_mean

def find_average_position_shifted(realwf_array, procent, center):
    """
    Finds the average position of the wavefuntion (weighted center of mass)
    Input:
        realwf_array: array with realspace wavefuntions as output by vaspwfc.wfc_r
        procent: sets psi=0 when |psi|^2 < max * procent
        center:
    Returns:
        average position
    """
    sh = realwf_array[0].shape
    wf2_tot = np.zeros(sh)

    for wf in realwf_array:
        wf2_tot += wf*wf.conj()

    cutoff = np.max(wf2_tot)*procent


    [a1,a2,a3] = np.meshgrid(np.arange(sh[0])/(sh[0]-1),np.arange(sh[1])/(sh[1]-1),np.arange(sh[2])/(sh[2]-1), indexing="ij")

    wf2_cut = np.where(wf2_tot > cutoff, wf2_tot, 0)
    wf2_cut = wf2_cut/np.sum(wf2_cut)

    c = np.array(center)-0.5
    def shift(arr, *c):
        for i, a in enumerate(arr):
            if a-c >= 0.5:
                arr[i] -= 1
            elif a-c < -0.5:
                arr[i] += 1
        return arr

    a1 = ((a1-c[0]) % 1 + c[0])
    a2 = ((a2-c[1]) % 1 + c[1])
    a3 = ((a3-c[2]) % 1 + c[2])


    avg_pos = [np.sum(a1*wf2_cut), np.sum(a2*wf2_cut), np.sum(a3*wf2_cut)]

    return np.array(avg_pos) % 1

def get_sg(sym, pos_file, settings):
    """
    Fetches space group and point group

    Input:
        sym: An object of the Symmetry class (from aflow_sym_python)
        pos_file: string that is the path to a crystal structure file like POSCAR or CONTCAR
        settings:

    Returns:
        Space group number
        Point group name
    """

    out = sym.get_sgdata(pos_file,settings['aflow_tolerance'])
    return out['space_group_number'], out['space_group_Schoenflies']



def get_energy_and_band_degen(eig_file,spin_i,lower_b,upper_b, settings):
    """
    Extracts band index, eigenvalue and occupation,
    and sorts by degeneracy

    Input:
        eig_file: path to EIGENVAL file
        spin_i: spin channel
        lower_b: index of lower band
        upper_b: index of upper band
        settings: path to settings file

    Returns:
        bands sorted by degeneracy
        band eigenvalues by degeneracy
        band occpation by degeneracy
    """

    # Get energies from EIGENVAL
    file = open(eig_file,"r")
    contents = file.readlines()
    file.close()

    band_num = []
    band_en = []
    band_occ = []

    for i in range(lower_b,upper_b+1):
        line = contents[i+7].split()
        band_num.append(line[0])
        if spin_i == 1:
            band_en.append(line[1])
            band_occ.append(float(line[3]))
        if spin_i == 2:
            band_en.append(line[2])
            band_occ.append(float(line[4]))

    bands_by_degen = []
    band_en_by_degen = []
    band_occ_by_degen = []
    temp = []
    occ = 0

    # Sort bands by pairing those with degenerate energies
    for i, E in enumerate(band_en):
        temp.append(int(band_num[i]))
        occ += band_occ[i]
        if i+1 >= len(band_en):
            bands_by_degen.append(temp)
            band_en_by_degen.append(float(E))
            band_occ_by_degen.append(int(occ))
            temp = []
            occ = 0
        elif abs(float(E) - float(band_en[i+1])) > settings['degeneracy_tolerance']:
            bands_by_degen.append(temp)
            band_en_by_degen.append(float(E))
            band_occ_by_degen.append(int(occ))
            temp = []
            occ = 0

    return bands_by_degen, band_en_by_degen, band_occ_by_degen #,band_num, band_en,

def group_name_conv(gname):
    """
    Converts string of space group name from AFLOW-SYM output
    to string that corresponds to .lis file names

    Input:
        gname: point group name from AFLOW-SYM

    Returns:
        point group name
    """

    l = list(gname)
    try:
        i = l.index('}')
        l = l[0:i]
    except:
        0
    try:
        l.remove('_')
        l.remove('{')
    except:
        0
    try:
        l.remove('\n')
    except:
        0

    str = ''
    for j in range(len(l)):
        str = str + l[j]
    if str == 'Ci':
        str = 'S2'
    if str == 'Cs' or str == 'CS':
        str = 'C1h'
    return str

def get_character_table(gname,settings):
    """
    Get character table as array from gname.lis

    Input:
        gname: name of point group
        settings: path to settings file

    Returns:
        character table
        position vectors and their representation
    """
    file = open(settings['char_table_dir']+"/"+gname+".lis","r")

    lines = file.readlines()
    endoftable = lines[1:-1].index(" \n")
    lines = lines[1:endoftable+1]

    topline_list = lines[0].split()
    for i, str in enumerate(topline_list):
        if list(str)[0] == '<':
            topline_list = topline_list[0:i]
            break

    last_col = []
    for row in lines[1:]:
        for i, str in enumerate(row.split()[1:]):
            if list(str)[0] == '.' or list(str)[0] == 'T':
                last_col.append(i+1)
                break

    char_table = [topline_list]

    for i in range(1,len(lines)):
        row = lines[i].split()[0:last_col[i-1]]
        try:
            row.remove("*")
            last_col[i-1] = last_col[i-1] - 1
            print(row[0]+" is reducible.")
        except:
            0
        char_table.append(row)

    pos_vec_rep = []
    for i, row in enumerate(lines[1:]):
        row = row.split()
        str = list(row[last_col[i]+1])
        for s in str:
            if s == "T":
                pos_vec_rep.append([char_table[i+1],str])
                break

    for i, component in enumerate(pos_vec_rep):
        comp = component[1]
        cartesian_comp_str = ""
        for n in range(3):
            cart_sym = ""
            if comp[n] == 'T':
                if n == 0:
                    cart_sym = " x"
                elif n == 1:
                    cart_sym = " y"
                elif n == 2:
                    cart_sym = " z"
            cartesian_comp_str = cartesian_comp_str + cart_sym
        pos_vec_rep[i][1] = cartesian_comp_str


    return char_table, pos_vec_rep


def get_vb_and_cb(eig_file, lower_b, upper_b, spin_channel):
    """
    Returns eigenvalues (energies) of the bands just below and above the
    bands considered in the symmetry analysis.

    Inputs:
        eig_file: path to EIGENVAL files
        lower_b: lowest band analysed (should be first above vb)
        upper_b: highest band analysed (should be first below cb)
        spin_channel: 1 or 2 for spin up or down (should make little differance)
    Returns:
        valance band eigenvalue
        conduction band eigenvalue
    """

    f = open(eig_file, "r")
    lines = f.readlines()
    f.close()

    if type(lower_b) == type([]):
        lower_b = lower_b[0]
    vb = lines[lower_b+6].split()[spin_channel]

    if type(upper_b) == type([]):
        upper_b = upper_b[-1]
    cb = lines[upper_b+8].split()[spin_channel]

    return float(vb), float(cb)

def get_transformation_matrix(pos_file):
    """
    Inputs:
        pos_file: string with path to POSCAR or CONTCAR containing lattice vectors
    Returns:
        3x3 transformation matrix that converts cartesian coordinates to lattice basis coordinates
    """

    f = open(pos_file, "r")
    lines = f.readlines()[2:5]
    f.close()

    vectors = []
    temp_row = []

    for i in range(3):
        for j in range(3):
            temp_row.append(float(lines[j].split()[i]))
        vectors.append(temp_row)
        temp_row = []

    return vectors

def latticel2norm(vector,pos_file):
    """
    Calculates the cartesian length of a vector given in the lattice basis
    Inputs:
        vector: numpy vector in lattice basis
        pos_file: string with path to POSCAR or CONTCAR containing lattice vectors
    Returns:
        proper L2-norm, length in Ångström
    """

    T = get_transformation_matrix(pos_file)
    cart_vec = T*vector

    return np.linalg.norm(cart_vec)

def gvectors_and_energy(wav, ikpt=1, force_Gamma=False, check_consistency=False):
    '''
    This is a modified version of the gvectors() function in vaspwfc.

    Generate the G-vectors that satisfies the following relation
        (G + k)**2 / 2 < ENCUT
    Also outputs the energy corresponding to those gvectors

    Inputs:

    '''
    assert 1 <= ikpt <= wav._nkpts,  'Invalid kpoint index!'
    kvec = wav._kvecs[ikpt-1]
    # force_Gamma: consider gamma-only case regardless of the actual setting
    lgam = True if force_Gamma else wav._lgam

    # fx, fy, fz = [fftfreq(n) * n for n in wav._ngrid]
    # fftfreq in scipy.fftpack is a little different with VASP frequencies
    ############################################################
    # Gamma version -50% memory usage and 1x speed.
    ############################################################
    # fx = [ii if ii < wav._ngrid[0] // 2 + 1 else ii - wav._ngrid[0]
    #       for ii in range(
    #           wav._ngrid[0] // 2 + 1
    #           if (lgam and (wav._gam_half == 'x'))
    #           else
    #           wav._ngrid[0])]
    # fy = [jj if jj < wav._ngrid[1] // 2 + 1 else jj - wav._ngrid[1]
    #       for jj in range(wav._ngrid[1])]
    # fz = [kk if kk < wav._ngrid[2] // 2 + 1 else kk - wav._ngrid[2]
    #       for kk in range(
    #           wav._ngrid[2] // 2 + 1
    #           if (lgam and (wav._gam_half == 'z'))
    #           else
    #           wav._ngrid[2])]

    fx, fy, fz = [np.arange(n, dtype=int) for n in wav._ngrid]
    fx[wav._ngrid[0] // 2 + 1:] -= wav._ngrid[0]
    fy[wav._ngrid[1] // 2 + 1:] -= wav._ngrid[1]
    fz[wav._ngrid[2] // 2 + 1:] -= wav._ngrid[2]
    if lgam:
        if wav._gam_half == 'x':
            fx = fx[:wav._ngrid[0] // 2 + 1]
        else:
            fz = fz[:wav._ngrid[2] // 2 + 1]

    # if lgam:
    #     # parallel gamma version of VASP WAVECAR exclude some planewave
    #     # components, -DwNGZHalf
    #     if wav._gam_half == 'z':
    #         kgrid = np.array([(fx[ii], fy[jj], fz[kk])
    #                           for kk in range(wav._ngrid[2])
    #                           for jj in range(wav._ngrid[1])
    #                           for ii in range(wav._ngrid[0])
    #                           if (
    #                               (fz[kk] > 0) or
    #                               (fz[kk] == 0 and fy[jj] > 0) or
    #                               (fz[kk] == 0 and fy[jj]
    #                                == 0 and fx[ii] >= 0)
    #         )], dtype=float)
    #     else:
    #         kgrid = np.array([(fx[ii], fy[jj], fz[kk])
    #                           for kk in range(wav._ngrid[2])
    #                           for jj in range(wav._ngrid[1])
    #                           for ii in range(wav._ngrid[0])
    #                           if (
    #                               (fx[ii] > 0) or
    #                               (fx[ii] == 0 and fy[jj] > 0) or
    #                               (fx[ii] == 0 and fy[jj]
    #                                == 0 and fz[kk] >= 0)
    #         )], dtype=float)
    # else:
    #     kgrid = np.array([(fx[ii], fy[jj], fz[kk])
    #                       for kk in range(wav._ngrid[2])
    #                       for jj in range(wav._ngrid[1])
    #                       for ii in range(wav._ngrid[0])], dtype=float)

    ############################################################
    # 10x faster
    ############################################################
    # In meshgrid, fx run the fastest, fz the slowest
    gz, gy, gx = np.array(
        np.meshgrid(fz, fy, fx, indexing='ij')
    ).reshape((3, -1))
    kgrid = np.array([gx, gy, gz], dtype=float).T
    if lgam:
        if wav._gam_half == 'z':
            kgrid = kgrid[
                (gz > 0) |
                ((gz == 0) & (gy > 0)) |
                ((gz == 0) & (gy == 0) & (gx >= 0))
            ]
        else:
            kgrid = kgrid[
                (gx > 0) |
                ((gx == 0) & (gy > 0)) |
                ((gx == 0) & (gy == 0) & (gz >= 0))
            ]

    # Kinetic_Energy = (G + k)**2 / 2
    # HSQDTM    =  hbar**2/(2*ELECTRON MASS)
    KENERGY = HSQDTM * np.linalg.norm(
        np.dot(kgrid + kvec[np.newaxis, :], TPI*wav._Bcell), axis=1
    )**2
    # find Gvectors where (G + k)**2 / 2 < ENCUT
    Gvec = kgrid[np.where(KENERGY < wav._encut)[0]]

    # Check if the calculated number of planewaves and the one recorded in the
    # WAVECAR are equal
    if check_consistency:
        # if wav._lsoc:
        #     assert Gvec.shape[0] == wav._nplws[ikpt - 1] // 2, \
        #         'No. of planewaves not consistent for an SOC WAVECAR! %d %d %d' % \
        #         (Gvec.shape[0], wav._nplws[ikpt - 1],
        #          np.prod(wav._ngrid))
        # else:
        #     assert Gvec.shape[0] == wav._nplws[ikpt - 1], 'No. of planewaves not consistent! %d %d %d' % \
        #         (Gvec.shape[0], wav._nplws[ikpt - 1],
        #          np.prod(wav._ngrid))

        if Gvec.shape[0] != wav._nplws[ikpt - 1]:
            if Gvec.shape[0] * 2 == wav._nplws[ikpt - 1]:
                if not wav._lsoc:
                    raise ValueError('''
                    It seems that you are reading a WAVECAR from a NONCOLLINEAR VASP.
                    Please set 'lsorbit = True' when loading the WAVECAR.
                    For example:

                        wfc = vaspwfc('WAVECAR', lsorbit=True)
                    ''')
            elif Gvec.shape[0] == 2 * wav._nplws[ikpt - 1] - 1:
                if not wav._lgam:
                    raise ValueError('''
                    It seems that you are reading a WAVECAR from a GAMMA-ONLY VASP.  Please set
                    'lgamma = True' when loading the WAVECAR.  Moreover, you may want to set
                    "gamma_half" if you are using VASP version <= 5.2.x.  For VASP <= 5.2.x, check
                    which FFT VASP uses by the following command:

                        $ grep 'use.* FFT for wave' OUTCAR

                    Then

                        # for parallel FFT, VASP <= 5.2.x
                        wfc = vaspwfc('WAVECAR', lgamma=True, gamma_half='z')

                        # for serial FFT, VASP <= 5.2.x
                        wfc = vaspwfc('WAVECAR', lgamma=True, gamma_half='x')

                    For VASP >= 5.4, WAVECAR is written with x-direction half grid regardless of
                    parallel or serial FFT.

                        # "gamma_half" default to "x" for VASP >= 5.4
                        wfc = vaspwfc('WAVECAR', lgamma=True, gamma_half='x')
                    ''')
            else:
                raise ValueError('''
                NO. OF PLANEWAVES NOT CONSISTENT:

                    THIS CODE -> %d
                    FROM VASP -> %d
                       NGRIDS -> %d
                ''' % (Gvec.shape[0],
                       wav._nplws[ikpt - 1] // 2 if wav._lsoc else wav._nplws[ikpt - 1],
                       np.prod(wav._ngrid))
                )

    return np.asarray(Gvec, dtype=int), KENERGY

def openf(filename):
    """
    Opens file in read mode, either regular or compressed
    Inputs:
        filename: name of file
    Returns:
        file
    """

    if filename[-3:] == "bz2":
        file = bz2.open(filename, "r")
    else:
        file = open(filename, "r")

    return file

def load_settings(file):
    """
    Settings file for both overlap calculations and symmetry analysis.
    Input:
        file: string with path to settings file
    Returns:
        dicitonary with settings
    """

    try:
        f = open(file,"r")
        data = json.load(f)
        f.close()
        print("Loaded: ",str(file))
    except:
        data = {}
        print("Loaded default settings.")

    try:
        data['aflow_tolerance']
    except Exception as e:
        data['aflow_tolerance'] = "tight"
    try:
        data['degeneracy_tolerance']
    except Exception as e:
        data['degeneracy_tolerance'] = 0.001
    try:
        data['round_if_close']
    except Exception as e:
        data['round_if_close'] = True
    try:
        data['round_if_close_tolerance']
    except Exception as e:
        data['round_if_close_tolerance'] = 0.05
    try:
        data['irrep_tolerance']
    except Exception as e:
        data['irrep_tolerance'] = 0.05
    try:
        data['tdm_irrep_from_irrep']
    except Exception as e:
        data['tdm_irrep_from_irrep'] = True
    try:
        data['Gvec_reduction']
    except Exception as e:
        data['Gvec_reduction'] = 0.2
    try:
        data['Gammapoint_calc']
    except Exception as e:
        data['Gammapoint_calc'] = True
    try:
        data['realgrid_mult']
    except Exception as e:
        data['realgrid_mult'] = 4
    try:
        data['percent_cutoff']
    except Exception as e:
        data['percent_cutoff'] = 0.40
    try:
        data['char_table_dir']
    except Exception as e:
        data['char_table_dir'] = "/dedur01/data/wilst89/software/character_tables"
    return data

if __name__ == "__main__":
    0
