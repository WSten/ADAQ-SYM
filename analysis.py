#!/usr/bin/env python3

from vaspwfc import vaspwfc
import numpy as np
import math as m
import os
from extract import *
from overlap import *
from pprint import pprint
import pickle

def get_character(eig_file, folder_path_out, name, settings, spin_i,lower_b,upper_b,permutation):
    """
    Inputs:
        eig_file: string that is the path to a EIGENVAL file
        folder_path_out: string that is the path to output directory
        name: string with name (numbering) of defect
        settings_file: string that is the path to a json file
        spin_i: 1 or 2 for spin up or down
        lower_b: index of lowest considered band
        upper_b: index of highest considered band
        permutation: permutation of which symmetry belongs to each class

    Returns:
        characters
        band index
        band eigenvalue
        band occupation, all output grouped by degeneracy
    """

    ov_path = os.path.join(folder_path_out,"Overlaps_"+name+".pickle")
    f = open(ov_path,"rb")
    sym_ov_array = pickle.load(f)
    f.close()

    round_if_close = settings['round_if_close']
    tol = settings['round_if_close_tolerance']

    bands_by_degen, band_en_by_degen, band_occ_by_degen = get_energy_and_band_degen(eig_file,spin_i,lower_b,upper_b,settings)
    overlap_list = []


    for deg_bands in bands_by_degen:
        temp = []
        for sym_i in range(len(sym_ov_array[0][1])):
            sum = 0

            for band_nr in deg_bands:
                band_i = int(band_nr)-lower_b
                sum = sum + sym_ov_array[band_i][1][sym_i][4]
            temp.append(sum)
        overlap_list.append(temp)
    #print("overlap summed")
    #pprint(overlap_list)
    #print("____")
    sym_list = []

    for sym_i in range(len(sym_ov_array[0][1])):
        sym_list.append(sym_ov_array[0][1][sym_i][1])

    sym_list, overlap_list = average_by_class(permutation, sym_list,overlap_list,bands_by_degen)
    #print("After truncate")
    #pprint(sym_list)
    #pprint(overlap_list)
    #print("___________")
    character_list = []

    for deg_band_i in range(len(overlap_list)):
        temp = []
        for sym_i in range(len(overlap_list[0])):
            ov = overlap_list[deg_band_i][sym_i]
            if abs(np.imag(ov)) > tol:
                print(deg_band_i,spin_i, sym_i)
                print(ov)
                print("Imaginary part of overlap > "+str(tol)+" !!!")
                print("_______")
                file = open("Imag_char.txt","a+")
                file.write(str(bands_by_degen[deg_band_i])+"  "+str(spin_i)+"  "+str(sym_ov_array[deg_band_i][1][sym_i][1])+"  "+str(ov)+"\n")
                file.close()
            ov_r = np.real(ov)
            if abs(round(ov_r)-ov_r) < tol:
                if round_if_close and abs(np.imag(ov)) < tol:
                    temp.append(round(ov_r))
                else:
                    temp.append(ov)
            else:
                temp.append(ov)
        character_list.append(temp)

    characters_by_degen = []
    characters_by_degen.append(sym_list)
    for i in range(len(character_list)):
        characters_by_degen.append(character_list[i])

    return characters_by_degen, bands_by_degen, band_en_by_degen, band_occ_by_degen


def average_by_class(permutation, symmetry_list,overlap_list,bands_by_degen):
    """
    Takes the average of the overlaps within each congucay class.

    Inputs:
        permutation: permutation of which symmetry belongs to each class
        symmetry_list: list of symbols for each symmetry operator
        overlap_list: list of overlaps for each symmetry operator
        bands_by_degen: bands grouped by degeneracy

    Returns:
        list with the symbol for each class
        list of average overlap for each class
    """
    sym_list_trunc = []
    for l in permutation:
        sym_list_trunc.append(symmetry_list[l[0]])

    overlap_list_trunc = []
    temp_ov = []
    sum = 0

    for deg_band_i in range(len(bands_by_degen)):
        for l in permutation:
            for i in l:
                sum += overlap_list[deg_band_i][i]
            temp_ov.append(sum/len(l))
            sum = 0
        overlap_list_trunc.append(temp_ov)
        temp_ov = []

    return sym_list_trunc, overlap_list_trunc

def order_columns_old(Sym_ops, name, settings, ch_table):


    #sym = Symmetry()
    #s = get_symmetry_operators(sym, host_pos_file, settings)
    #T = get_transformation_matrix(pos_file)

    sym_ops = Sym_ops[1]

    #for Sym in s[1]:
    #    U = np.linalg.inv(T).dot(np.array(Sym).dot(T))
    #    U = round_matrix(U,settings)
    #    sym_ops.append(U)


    class_symbols, perm, mult = sort_into_classes(sym_ops, Sym_ops[0])
    #print(s[0])

    symbol_table = ch_table[0][1:]

    for m in mult:
        if m != 1:
            symbol_table.remove(str(m))
    print(class_symbols)
    print(symbol_table)
    c_perm = []

    for symb in symbol_table:
        symb_2 = ""
        if symb == "E":
            symb = "1"
        elif symb == "sv" or symb == "sd" or symb == "sh":
            symb = "S2"
            symb_2 = "s"
        elif symb == "S6": # !!!!???? S_n vs S_2n
            symb = "S3"
        elif symb == "C2'" or symb == "C2''" or symb == 'C2"':
            symb = "C2"
        print(symb, symb_2)
        try:
            c_perm.append(class_symbols.index(symb))

        except:
            c_perm.append(class_symbols.index(symb_2))

        print(class_symbols)

    print(class_symbols)
    print(c_perm)
    permutation = []
    multiplicity = []

    for i in c_perm:
        permutation.append(perm[i])
        multiplicity.append(len(perm[i]))
    print(permutation)
    print(multiplicity)
    #return sym_ops, permutation, multiplicity
    return permutation, multiplicity

def order_columns(Sym_ops, name, settings, ch_table):
    """
    Orders the columns (lists) according to how character tables are ordered.

    Inputs:
        Sym_ops: array with symmetry operator info,
                 output for get_symmetry_operators()
        name: string with name (numbering) of defect
        settings: settings dicitonary
        ch_table: character table as array

    Returns:
        permutation of which symmetry belongs to each class
        multiplicity (order) of each class
    """

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

    class_symbols, class_perm, mult = sort_into_classes(Sym_ops[1], Sym_ops[0])

    perm = []

    for i, symb in enumerate(classes_right_order):
        symb_2 = "X"
        if symb == "E":
            symb = "1"
        elif symb == "sv" or symb == "sd" or symb == "sh":
            symb = "S2"
            symb_2 = "s"
        elif symb == "S6": # !!!!???? S_n vs S_2n
            symb = "S3"
        elif symb == "C2'" or symb == "C2''" or symb == 'C2"':
            symb = "C2"
        occurance = [j for j, n in enumerate(class_symbols) if n == symb or n == symb_2]
        for ind in occurance:
            if mult_right_order[i] == mult[ind] and not ind in perm:
                perm.append(ind)

    permutation = []
    for p in perm:
        permutation.append(class_perm[p])

    return permutation, mult_right_order


def get_irrep_symbols(ch_table):
    """
    Get irrep symbols from first column of character table
    """
    symb = []
    for row in ch_table[1:]:
        symb.append(row[0])
    return symb


def get_rep(ch_table, chars, mult, irrep_tol):
    """
    Get representation from characters by projecting the characters on the
    characters of the irreps.

    Inputs:
        ch_table: character table as Array
        chars: array of characters
        mult: multiplicity (order) of each class
        irrep_tol: tolerance
    Returns:
        representation as array
    """

    # Through projection of irreps
    irreps = []
    #ch = np.asfarray(chars,float)
    ch = chars

    for row in ch_table[1:]:
        irr = row[0]
        row = np.array(list(int(r) for r in row[1:]))
        proj = np.dot(mult*row,ch)/sum(mult)
        #print(dot_p)

        proj_r = np.real(proj)
        if abs(np.imag(proj)) > irrep_tol:
            irreps.append(0)
        elif abs(proj_r-round(proj_r)) < irrep_tol:
            irreps.append(round(proj_r))
        else:
            irreps.append(0)
    if irreps == []:
        irreps = np.zeros(len(ch_table[1:]))
    return irreps


def get_total_irrep_string(irr_symbols, irreps):
    """
    Get string with total representation
    Inputs:
        irr_symbols: list of irrep symbols as strings
        irreps: array with irreps
    Returns:
        representation string
    """
    string = ""
    count = 0
    for i, int in enumerate(irreps):
        if int == 1:
            if count != 0:
                string = string + "+ "
            string = string + irr_symbols[i]+" "
            count += 1
        elif int != 0:
            if count != 0:
                string = string + "+ "
            count += 1
            string = string + str(int) + irr_symbols[i]+" "

    return string

def get_csm(ch_table, chars, mult, irrep_tol):
    """
    Get contious symmetry measure:
    How close a set of characters is to an irrep
    0 = Symmetric, 100 = Asymmetric
    111.1 = Error, imaginary component too large
    Inputs:
        ch_table: character table as Array
        chars: array of characters
        mult: multiplicity (order) of each class
        irrep_tol: tolerance
    Returns:
        representation as array
    """

    irrep_csm = []
    ch = chars

    for row in ch_table[1:]:
        irr = row[0]
        row = np.array(list(int(r) for r in row[1:]))
        proj = np.dot(mult*row,ch)/sum(mult)

        proj_r = np.real(proj)
        if abs(np.imag(proj)) > irrep_tol:
            irrep_csm.append(111.1)
        else:
            irrep_csm.append(round(100-100*proj_r,1))
    if irrep_csm == []:
        irrep_csm = 100*np.ones(len(ch_table[1:]))
    return irrep_csm

def gather_csm(eig_file, name, spin_i, lower_b, upper_b, Sym_ops, PGname, folder_path_out, settings):
    """
    Finds csm for considered bands
    Input:
        eig_file: path to EIGENVAL file
        name: string with name
        spin_i: 1 or 2 for spin up or down
        lower_b: index of lowest considered band
        upper_b: index of highest considered band
        Sym_ops: symmetry operations in array
        PGname: string with name of point group
        folder_path_out: path of output directory
        settings: settings dicitonary
    Returns:
        csm_array
        irreps
    """

    irrep_tol = settings['irrep_tolerance']

    ch_table, pos_vector = get_character_table(PGname,settings)

    perm, mult = order_columns(Sym_ops, name, settings, ch_table)

    ch_list, deg, en_deg, occ_deg = get_character(eig_file,folder_path_out,name,settings, spin_i,lower_b,upper_b,perm)

    irreps = [ch[0] for ch in ch_table[1:]]
    irreps = get_irrep_symbols(ch_table)
    csm_array = []

    for i in range(len(deg)):
        csm_array.append([deg[i], get_csm(ch_table, ch_list[1+i], mult, irrep_tol)])

    return csm_array, irreps

def csm_main(s1bands, s2bands, PGname, Sym_ops, settings, eig_file = "EIGENVAL",  name="", folder_path_out=""):
    """
    Calculates csm and writes to a file.
    Input:
        s1bands: list of band indices for spin channel 1
        s2bands: list of band indices for spin channel 2
        PGname: name of point group
        Sym_ops: array with symmetry operator info,
                 output for get_symmetry_operators()
        settings: settings dicitonary
        eig_file: string that is the path to a EIGENVAL file
        name: string with name (numbering) of defect
        folder_path_out: string that is the path to output directory
    Returns:

    """

    csm_array_s1 = []
    lower_b_s1 = s1bands[0]
    upper_b_s1 = s1bands[-1]
    name1 = name+"_S1"

    csm_array_s2 = []
    lower_b_s2 = s2bands[0]
    upper_b_s2 = s2bands[-1]
    name2 = name+"_S2"

    if len(s1bands) > 0:
        csm_array_s1, irreps = gather_csm(eig_file, name1, 1, lower_b_s1, upper_b_s1, Sym_ops, PGname, folder_path_out, settings)

    if len(s2bands) > 0:
        csm_array_s2, irreps = gather_csm(eig_file, name2, 2, lower_b_s2, upper_b_s2, Sym_ops, PGname, folder_path_out, settings)

    csm_path = os.path.join(folder_path_out,"CSM"+name+".txt")
    file = open(csm_path,"w+")
    file.write("Continous Symmetry Measure\n")
    file.write("0 = Symmetric, 100 = Fully Asymmetric\n")
    file.write("111.1 = Error, imaginary component too large\n")
    file.write("\nPoint group: "+PGname+"\n")

    file.write("Spin Up: \n")
    irrep_string = f"{'Band:':<20}"
    for irr in irreps:
        irrep_string += f"{irr:<7}"
    irrep_string += "\n"
    file.write(irrep_string)


    for csm in csm_array_s1:
        csm_string = f"{str(csm[0]):<20}"
        for c in csm[1]:
            csm_string += f"{c:<7}"
        file.write(csm_string+"\n")
    file.write("\n\nSpin Down: \n")
    file.write(irrep_string)
    for csm in csm_array_s2:
        csm_string = f"{str(csm[0]):<20}"
        for c in csm[1]:
            csm_string += f"{c:<7}"
        file.write(csm_string+"\n")

    file.close()
    return 0

def get_allowed_transitions(char_table, position_vecs, char_list, multiplicity, bands_by_degen, eigen_by_bands, occ_by_degen, settings):
    """
    Takes table of characters and calculates the representation of each possible
    transition for each polarisation. Checks if representation contains
    identity representation.

    Input:
        char_table: character table as numpy array
        position_vecs: characters of linear functions read from character tables
        char_list: array with characters of each bands
        multiplicity: order of each conjugacy class of the point group_name_conv
        bands_by_degen: bands ordered by degeneracy
        eigen_by_bands: eigenvalues of each band, according to above order
        occ_by_degen: occupation of each band, according to above order
        settings: settings dicitonary
    Returns:
        Array with information on each transition
        list of bands where no irrep was found
    """

    irrep_tol = settings['irrep_tolerance']
    tdm_irrep_from_irrep = settings['tdm_irrep_from_irrep']
    ch_list = char_list[1:]
    irr_symbols = get_irrep_symbols(char_table)
    irreps = []
    no_irr = []

    for i, ch in enumerate(ch_list):
        irr = get_rep(char_table,ch,multiplicity,irrep_tol)
        irreps.append(irr)
        irr = np.array(irr)
        if np.sum(irr*irr) == 0:
            for band in bands_by_degen[i]:
                no_irr.append(band)

        if tdm_irrep_from_irrep:
            ch_list[i] = np.array([round(np.real(c)) for c in ch])

    transitions = []

    for i, ch_init in enumerate(ch_list):

        for j, ch_final in enumerate(ch_list):

            #if i < j:
            if i != j and 0 < occ_by_degen[i] and occ_by_degen[j] < len(bands_by_degen[j]):
                for pos_vec in position_vecs:
                    pos_v = np.asfarray(pos_vec[0][1:],float)

                    tdm = ch_final * (pos_v * ch_init)

                    tdm_irr = get_rep(char_table, tdm, multiplicity,irrep_tol)

                    # For saving with pickle
                    tdm = tdm.tolist()

                    bool = (tdm_irr[0] >= 1)

                    transitions.append([bands_by_degen[i], eigen_by_bands[i], get_total_irrep_string(irr_symbols,irreps[i]), \
                    bands_by_degen[j], eigen_by_bands[j], get_total_irrep_string(irr_symbols,irreps[j]), \
                    pos_vec[1], pos_vec[0][0], tdm, get_total_irrep_string(irr_symbols,tdm_irr), \
                    bool])

    return transitions, no_irr

def analyse_symmetry(eig_file, name, spin_i, lower_b, upper_b, Sym_ops, PGname, folder_path_out, settings):
    """
    Finds characters and irreps of the considered bands,
    finds allowed transitions and writes this to a file,
    also saves an array of info as a .npy file.
    Input:
        eig_file: path to EIGENVAL file
        name: string with name
        spin_i: 1 or 2 for spin up or down
        lower_b: index of lowest considered band
        upper_b: index of highest considered band
        Sym_ops: array with symmetry operator info,
                 output for get_symmetry_operators()
        PGname: string with name of point group
        folder_path_out: path of output directory
        settings: settings dicitonary
    Returns:
        list of bands where no irrep was found
    """

    irrep_tol = settings['irrep_tolerance']

    ch_table, pos_vector = get_character_table(PGname,settings)

    perm, mult = order_columns(Sym_ops, name, settings, ch_table)

    ch_list, deg, en_deg, occ_deg = get_character(eig_file,folder_path_out,name,settings, spin_i,lower_b,upper_b,perm)

    print(ch_list)
    pprint(deg)
    print("\n")

    tr_path = os.path.join(folder_path_out,"Transitions_"+name+".txt")
    file = open(tr_path,"w+")

    file.write("Character table:") #+str(ch_table[0])+str(ch_table[1])+str(ch_table[2])+str(ch_table[3])+"\n")
    for i in range(len(ch_table)):
        file.write("\n"+str(ch_table[i]))
    file.write("\nPosition vector representation:")# \n"+str(pos_vector[0])+"\n"+str(pos_vector[1]))
    for i in range(len(pos_vector)):
        file.write("\n"+str(pos_vector[i]))

    file.write(f"\n{'Band':<15} {'Eig':<10} {'Occ':<5} {'Rep':<5} {'Characters'} \n")

    band_info = []

    for i in range(len(deg)):

        tot_irr_string = get_total_irrep_string(get_irrep_symbols(ch_table),get_rep(ch_table, ch_list[1+i], mult,irrep_tol))
        print(i, deg[i], ch_list[1+i], tot_irr_string)
        #file.write(str(deg[i])+"   "+str(en_deg[i])+"   "+str(occ_deg[i])+"   "+str(ch_list[1+i])+"   "+tot_irr_string+"\n")
        file.write(f"{str(deg[i]):<15} {en_deg[i]:<10} {occ_deg[i]:<5} {tot_irr_string:<5} {ch_list[1+i]}\n")
        band_info.append([deg[i], en_deg[i], occ_deg[i], ch_list[1+i], tot_irr_string])



    tr, no_irr = get_allowed_transitions(ch_table, pos_vector, ch_list, mult, deg, en_deg, occ_deg, settings)
    #pprint(tr)

    #tr_array = np.asanyarray([PGname, band_info, tr],dtype=object)
    tr_array = [PGname, band_info, tr]
    tr2_path = os.path.join(folder_path_out,"Transitions_"+name+".pickle")
    #np.save(tr2_path, tr_array)
    f = open(tr2_path,"wb")
    pickle.dump(tr_array,f,protocol=2)
    f.close()


    file.write("\nInitial band, eigenvalue and irrep   Final band, eigenvalue and irrep    Position vector components and irrep  \n  TDM character and rep    Transition allowed?")

    for line in tr:
        file.write("\n"+str(line[0])+"  "+str(line[1])+"  "+str(line[2])+"  "+str(line[3])+\
        "  "+str(line[4])+"  "+str(line[5])+"  "+str(line[6])+"  "+str(line[7]))
        file.write("\n"+str(line[8])+"  "+str(line[9])+"  "+str(line[10])+"\n")

    file.close()



    return no_irr


def main(s1bands, s2bands, pos_file = "CONTCAR", wf_file = "WAVECAR", eig_file = "EIGENVAL", name="", settings_file = "settings.json", folder_path_out=""):
    """

    Input:
        s1bands: list of band indices for spin channel 1
        s2bands: list of band indices for spin channel 2
        pos_file: string that is the path to a crystal structure file like POSCAR or CONTCAR
        wf_file: string that is the path to a WAVECAR file
        eig_file: string that is the path to a EIGENVAL file
        name: string with name (numbering) of defect
        settings_file: string that is the path to a json file
        folder_path_out: string that is the path to output directory
    Returns:

    """


    settings = load_settings(settings_file)

    PGname = get_pointgroup(pos_file, settings)
    print("Point group: ",PGname)
    sym = Symmetry()
    Sym_ops = get_symmetry_operators(sym, pos_file, PGname, settings)

    """
    #Use this to manually reduce point group.
    temp = []
    for i in range(4):
        temp.append([Sym_ops[i][0], Sym_ops[i][3], Sym_ops[i][7], Sym_ops[i][11]])
    Sym_ops = temp
    PGname = "C2h"
    """

    # Initial overlap and analysis for spin up channel
    no_irr_s1 = []
    if len(s1bands) > 0:
        name1 = name+"_S1"
        lower_b_s1 = s1bands[0]
        upper_b_s1 = s1bands[-1]
        bands_by_degen_s1, band_en_by_degen_s1, band_occ_by_degen_s1 = get_energy_and_band_degen(eig_file,1,lower_b_s1,upper_b_s1,settings)
        print("Spin up orbitals by degeneracy: ", bands_by_degen_s1)
        try:
            center_path = os.path.join(folder_path_out,"Centers_"+name1+".npy")
            centers_s1 = np.load(center_path)
            print("Loaded centers!")
        except Exception as e:
            centers_s1 = get_orbital_centers(wf_file, bands_by_degen_s1, name1, 1, lower_b_s1, upper_b_s1, folder_path_out, settings)
            print("Calculated new centers!")
        print("Spin 1")
        res1 = get_overlaps_of_bands(wf_file, name1, 1, lower_b_s1, upper_b_s1, centers_s1, Sym_ops, folder_path_out, settings)
        write_overlaps_to_text(res1, folder_path_out, name1)
        no_irr_s1 = analyse_symmetry(eig_file, name1, 1, lower_b_s1, upper_b_s1, Sym_ops, PGname, folder_path_out, settings)
        good_centers_s1 = get_good_centers(name1, lower_b_s1, no_irr_s1, folder_path_out)

    # Initial overlap and analysis for spin down channel
    no_irr_s2 = []
    if len(s2bands) > 0:
        name2 = name+"_S2"
        lower_b_s2 = s2bands[0]
        upper_b_s2 = s2bands[-1]
        bands_by_degen_s2, band_en_by_degen_s2, band_occ_by_degen_s2 = get_energy_and_band_degen(eig_file,2,lower_b_s2,upper_b_s2,settings)
        print("Spin down orbitals by degeneracy: ", bands_by_degen_s2)
        try:
            center_path = os.path.join(folder_path_out,"Centers_"+name2+".npy")
            centers_s2 = np.load(center_path)
            print("Loaded centers!")
        except Exception as e:
            centers_s2 = get_orbital_centers(wf_file, bands_by_degen_s2, name2, 2, lower_b_s2, upper_b_s2, folder_path_out, settings)
            print("Calculated new centers!")
        print("Spin 2")
        res2 = get_overlaps_of_bands(wf_file, name2, 2, lower_b_s2, upper_b_s2, centers_s2, Sym_ops, folder_path_out, settings)
        write_overlaps_to_text(res2, folder_path_out, name2)
        no_irr_s2 = analyse_symmetry(eig_file, name2, 2, lower_b_s2, upper_b_s2, Sym_ops, PGname, folder_path_out, settings)
        good_centers_s2 = get_good_centers(name2, lower_b_s2, no_irr_s2, folder_path_out)



    # Gather good centers, add additional center candidates, e.g. impurity position
    good_centers = np.array([])
    if len(s1bands) > 0:
        good_centers = np.append(good_centers, good_centers_s1, axis=0)
    if len(s2bands) > 0:
        good_centers = np.append(good_centers, good_centers_s2, axis=0)

    defect_pos = get_single_species(pos_file)
    if defect_pos != None:
        good_centers = np.array([defect_pos])
    good_centers_s1 = good_centers
    good_centers_s2 = good_centers


    # Redo overlap and analysis with new centers if no irrep was found

    while 0 < len(no_irr_s1) and 0 < len(good_centers_s1):
        good_centers_s1, centers_s1 = replace_bad_centers(name1, lower_b_s1, no_irr_s1, good_centers_s1, folder_path_out)
        res1 = get_overlaps_of_bands(wf_file, name1, 1, lower_b_s1, upper_b_s1, centers_s1, Sym_ops, folder_path_out, settings)
        write_overlaps_to_text(res1, folder_path_out, name1)
        no_irr_s1 = analyse_symmetry(eig_file, name1, 1, lower_b_s1, upper_b_s1, Sym_ops, PGname, folder_path_out, settings)
    if len(no_irr_s1) > 0:
        file = open("no_irr.txt","a+")
        file.write("Spin1: "+str(no_irr_s1)+"\n")
        file.close()
    # Redo overlap and analysis with new centers if no irrep was found
    #if  0 < len(no_irr_s2) < len(s2bands):

    while 0 < len(no_irr_s2) and 0 < len(good_centers_s2):
        good_centers_s2, centers_s2 = replace_bad_centers(name2, lower_b_s2, no_irr_s2, good_centers_s2, folder_path_out)
        res2 = get_overlaps_of_bands(wf_file, name2, 2, lower_b_s2, upper_b_s2, centers_s2, Sym_ops, folder_path_out, settings)
        write_overlaps_to_text(res2, folder_path_out, name2)
        no_irr_s2 = analyse_symmetry(eig_file, name2, 2, lower_b_s2, upper_b_s2, Sym_ops, PGname, folder_path_out, settings)
    if len(no_irr_s2) > 0:
        file = open("no_irr.txt","a+")
        file.write("Spin2: "+str(no_irr_s2)+"\n")
        file.close()

    csm_main(s1bands, s2bands, PGname, Sym_ops, settings)

    return 0



if __name__ == "__main__":
    #host_pos_file="SiC/workflow_test/POSCAR_host"
    #main(1147,1149,2, host_pos_file="../hh/CONTCAR" ,name="test_6")

    main([], [1025, 1026, 1027], name="")
    #main([1153],[])

    print("Done")
