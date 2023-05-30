#!/usr/bin/env python

import os
import numpy as np
from aflow_sym_python import Symmetry
from extract import *
from overlap import *
from analysis import *
import pickle
import time
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import gridspec
import math as m
import sys

def fancy_subscript(string, lowercase):
    """
    Convert strings so matplotlib writes substripts.

    Inputs:
        string: irrep or point group
        lowercase: boolean to convert all letters to lowercase
    Returns:
        new string
    """
    if lowercase:
        string = string.lower()

    new_string = ''
    for s in string.split():
        new_string += s[0]
        if len(s) == 1:
            continue
        elif s[1] == "'":
            new_string += "\'"
        elif s[1] == '"':
            new_string += '"'
        else:
            new_string += '$_{'+s[1:]+'}$'

        new_string += ' '
    new_string  += ''
    return new_string

def plot_transitions(allowed_tr, spin, labels):
    """

    """

    # Set starting x-position of transition arrows
    if len(allowed_tr) == 1 and spin == 1:
        x_pos = 2.3
    elif spin == 1:
        x_pos = 0.5
    elif len(allowed_tr) == 1 and spin == 2:
        x_pos = 7.7
    else:
        x_pos = 6

    for i, tr in enumerate(allowed_tr):
        # Determine color and labels based on polarization
        if tr[6] == " x":
            c = 'r'
            p = '$\perp$'
        elif tr[6] == " y":
            c = 'g'
            p = '$\perp$'
        elif tr[6] == " z":
            c = 'b'
            p = '$\parallel$'
        elif tr[6] == " x y":
            c = 'r'
            p = '$\perp$'
        elif tr[6] == " x y z":
            c = 'r'
            p = ''

        labels.append((c,p+" "+fancy_subscript(tr[7], True)))

        # Draw arrow between initial and final band
        diff = float(tr[4])-float(tr[1])
        plt.arrow(x_pos,float(tr[1]), 0, 0.99*diff, color=c,\
                    length_includes_head=True, width=0.1, head_length = abs(diff)/10)

        if len(allowed_tr) != 1:
            x_pos += 3.5/(len(allowed_tr)-1)

    return labels

def generate_handles(labels):
    """

    """

    handles = []
    labels = list(set(labels))
    for c, p in labels:
        label_line = plt.vlines([], [], 0, color=c, label=p)
        handles.append(label_line)

    return handles

def plot_eigen(band_info, spin):
    """

    """

    # The following loops draw the eneregy levels and their occupation,
    # and also draws the irrep for each band.
    # One might need to manually adjust some of these positions to have a
    # readable diagram.

    for i, state in enumerate(band_info):
        eig = float(state[1])

        # Left and right bounds of energy level line
        if spin == 1:
            x_left = 0
            x_right = 4.5
        else:
            x_left = 5.5
            x_right = 10
        # Irrep position relative to x_right
        x_irrep = 0.1
        # Vertical offset
        y = 0.00

        if state[2] == 0:
            linestyle = "dotted"
        else:
            linestyle = "solid"

        # Draw energy level line
        plt.hlines(eig,x_left,x_right,linestyles=linestyle)
        width=0.06

        # Shift position of arrows indicating an occupied band
        shift = 0
        if i == 0:
            shift = 0
        if i == 1:
            shift = 0

        if state[2] == 2:
            shift = -0.5
        elif state[2] == 3:
            shift = -1

        # Print occupation arrows
        if spin == 1:
            for j in range(state[2]):
                plt.arrow(2.3+shift+j, eig-0.12, 0, 0.24, color='k', length_includes_head=False, \
                        width=width, head_width=3*width, head_length=0.1)
        else:
            for j in range(state[2]):
                plt.arrow(7.7+shift+j, eig+0.12, 0, -0.24, color='k', length_includes_head=False, \
                        width=width, head_width=3*width, head_length=0.1)

        # Manually adjust vertical position of irrep text
        y += 0.0
        if i == 0:
            y += 0.0
        if i == 1:
            y += -0.06

        # Print irrep
        plt.text(x_right+x_irrep,eig+y,fancy_subscript(str(state[4]), True),fontsize=15)

    return 0

def plot_ipr(HOB, eig_file, wf_file, gs):
    """

    """
    
    ax1 = plt.subplot(gs[0])
    plt.margins(0,0,tight=True)
    #plt.ylim(8.3, 14.3)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(True)
    ax1.spines['bottom'].set_visible(True)
    ax1.spines['left'].set_visible(False)
    ax1.invert_xaxis()
    ax1.yaxis.tick_right()
    #ax1.set_yticklabels([])

    wav = vaspwfc(wf_file, lgamma=True)

    # Grid density can be increased but this makes this script slower
    grid = wav._ngrid.copy() * 1

    # How many bands above and below HOB are considered for IPR
    extent = 15
    ipr_max = 0
    bands = range(HOB-extent,HOB+extent)
    for spin_channel in (1,2):
        iprs = []
        for b in bands:
            ks = [spin_channel,1,b]
            realwf = wav.wfc_r(*ks, ngrid=grid)
            wf_abs = np.abs(realwf)
            temp = np.sum(wf_abs**4)/np.sum(wf_abs**2)**2
            iprs.append(1000*temp) # IPR scaled by 1000 for readability

        ipr_max = np.max([ipr_max, np.max(iprs)])
        ipr_avg = np.average(iprs)

        print("Average IPR: ", ipr_avg)

        for i, ipr in enumerate(iprs):
            if ipr > ipr_avg:
                print("Band "+str(bands[i])+" Spin "+str(spin_channel)+" is localized. IPR: "+str(ipr))

        file = open(eig_file,"r")
        contents = file.readlines()
        file.close()
        eigenvalues = []
        for b in bands:
            eigenvalues.append(float(contents[b+7].split()[spin_channel]))

        if spin_channel == 1:
            c = 'm'
            spin_label = "Spin Up"
            avg_label = "Spin Up Avg"
        else:
            c = 'b'
            spin_label = "Spin Down"
            avg_label = "Spin Down Avg"
        plt.plot(iprs, eigenvalues, "x", color=c, label=spin_label)

        # Draw average IPR level in plot
        #plt.vlines(ipr_avg, eigenvalues[0], eigenvalues[-1], color=c, label=avg_label, linestyle='dotted')
        #plt.text(bands[-1]-5, ipr_avg*1.2, "Average IPR")

    #ax.set_ylim([9.05, 14.30])
    plt.xlabel("IPR [10^-3]", fontsize=15)
    plt.xlim([1.1*ipr_max, 0])
    #plt.legend(fontsize=15, bbox_to_anchor=(0.5, 0.90))
    plt.legend(fontsize=15)

    return ax1

def plot_levels_and_ipr(folder_path, plotname, HOB, eig_file= "EIGENVAL", filename="", wf_file="WAVECAR"):

    plt.rcParams['axes.linewidth'] = 2
    plt.rcParams['xtick.major.size'] = 5
    plt.rcParams['xtick.major.width'] = 1.5
    plt.rcParams['xtick.minor.size'] = 2.5
    plt.rcParams['xtick.minor.width'] = 1.5
    plt.rcParams['ytick.major.size'] = 5
    plt.rcParams['ytick.major.width'] = 1.5
    plt.rcParams['ytick.minor.size'] = 2.5
    plt.rcParams['ytick.minor.width'] = 1.5
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12
    plt.rcParams["figure.figsize"] = (14.4,7.2)

    HOB = int(HOB)

    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 3])
    plt.margins(0,0,tight=True)

    ax1 = plot_ipr(HOB, eig_file, wf_file, gs)

    plt.subplots_adjust(wspace=0.22)
    #ax2 = plt.subplot(1,2,2)
    ax2 = plt.subplot(gs[1], sharey=ax1)
    plt.margins(0,0,tight=True)

    plot_levels(folder_path, plotname)

    return 0

def plot_levels_one_spin(folder_path, filename, plotname, eig_file):

    plt.rcParams['axes.linewidth'] = 2
    plt.rcParams['xtick.major.size'] = 5
    plt.rcParams['xtick.major.width'] = 1.5
    plt.rcParams['xtick.minor.size'] = 2.5
    plt.rcParams['xtick.minor.width'] = 1.5
    plt.rcParams['ytick.major.size'] = 5
    plt.rcParams['ytick.major.width'] = 1.5
    plt.rcParams['ytick.minor.size'] = 2.5
    plt.rcParams['ytick.minor.width'] = 1.5
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12
    #plt.rcParams["figure.figsize"] = (8,4.8)

    #band_info, transition_info = np.load(folder_path+"/Transitions_"+filename+".npy",allow_pickle=True)
    tr_path = os.path.join(folder_path,"Transitions_"+filename+".pickle")
    f = open(tr_path, "rb")
    Pointgroup, band_info, transition_info = pickle.load(f)
    f.close()

    spin = int(list(filename)[-1])

    vb, cb = get_vb_and_cb(os.path.join(folder_path, eig_file), band_info[0][0], band_info[-1][0], spin)

    #plt.hlines(vb,0,10, linestyles='dashed')
    plt.text(-0.9,vb-0.22,"VB",fontsize=12)
    #plt.hlines(cb,0,10, linestyles='dashed')
    plt.text(-0.9,cb+0.02,"CB",fontsize=12)
    ax = plt.gca()
    vbrect = patches.Rectangle((0,vb),10,-0.6, color="0.7")
    cbrect = patches.Rectangle((0,cb),10,0.5, color="0.7")
    ax.add_patch(vbrect)
    ax.add_patch(cbrect)

    allowed_tr = []
    for tr in transition_info:
        if tr[-1]:
            allowed_tr.append(tr)



    ov_path = os.path.join(folder_path,"Overlaps_"+filename+".pickle")
    f = open(ov_path, "rb")
    ov = pickle.load(f)
    f.close()
    if len(ov[0][1]) == 1:
        principal_axis = np.array([0,0,0])
    else:
        principal_axis = np.array(ov[0][1][1][2])
        principal_axis = [round(p / max(principal_axis, key=abs), 3) for p in principal_axis]

    #vb = vb -0.75
    plt.text(0,cb+1.00, "Point group: "+fancy_subscript(Pointgroup, False), fontsize=12)
    plt.text(0,cb+0.65, "Principal axis: "+str(principal_axis), fontsize=12)

    labels = []
    labels = plot_transitions(allowed_tr, spin, labels)
    handles = generate_handles(labels)
    plot_eigen(band_info, spin)

    plt.xlim([-1, 15])
    plt.ylabel("Eigenvalue [eV]", fontsize=15)

    plt.xticks([])
    plt.box(on=None)
    plt.subplots_adjust(left=0.160)
    #plt.legend(handles = handles, bbox_to_anchor=(0.95, 0.04))
    plt.legend(handles = handles, bbox_to_anchor=(0.28, 0.72))

    plt.savefig(folder_path+"/Tr_"+plotname+".png")
    #plt.show()

    return 0

def plot_levels(folder_path, plotname, eig_file= "EIGENVAL", filename=""):

    plt.rcParams['axes.linewidth'] = 2
    plt.rcParams['xtick.major.size'] = 5
    plt.rcParams['xtick.major.width'] = 1.5
    plt.rcParams['xtick.minor.size'] = 2.5
    plt.rcParams['xtick.minor.width'] = 1.5
    plt.rcParams['ytick.major.size'] = 5
    plt.rcParams['ytick.major.width'] = 1.5
    plt.rcParams['ytick.minor.size'] = 2.5
    plt.rcParams['ytick.minor.width'] = 1.5
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12
    plt.rcParams["figure.figsize"] = (9.6,7.2)

    # Load point group, band and transition information
    tr_path = os.path.join(folder_path,"Transitions_"+filename+"_S1.pickle")
    f = open(tr_path, "rb")
    Pointgroup, band_info1, transition_info1 = pickle.load(f)
    f.close()
    tr_path = os.path.join(folder_path,"Transitions_"+filename+"_S2.pickle")
    f = open(tr_path, "rb")
    Pointgroup, band_info2, transition_info2 = pickle.load(f)
    f.close()

    # Get and draw vb and cb
    vb, cb = get_vb_and_cb(os.path.join(folder_path, eig_file), band_info2[0][0], band_info2[-1][0], 2)

    # Manual VB and CB for Diamond
    #vb=9.048191
    #cb=14.342933

    plt.text(-0.75,vb-0.35,"VB",fontsize=15)
    plt.text(-0.75,cb+0.15,"CB",fontsize=15)

    ax = plt.gca()
    vbrect = patches.Rectangle((0,vb),10,-0.7, color="0.7")
    cbrect = patches.Rectangle((0,cb),10,0.7, color="0.7")
    ax.add_patch(vbrect)
    ax.add_patch(cbrect)

    # Get allowed transitions
    allowed_tr1 = []
    allowed_tr2 = []
    for tr in transition_info1:
        if tr[-1]:
            allowed_tr1.append(tr)
    for tr in transition_info2:
        if tr[-1]:
            allowed_tr2.append(tr)

    # Read principal axis from overlap file
    ov_path = os.path.join(folder_path,"Overlaps_"+filename+"_S2.pickle")
    f = open(ov_path, "rb")
    ov = pickle.load(f)
    f.close()
    if len(ov[0][1]) == 1:
        principal_axis = np.array([0,0,0])
    else:
        principal_axis = np.array(ov[0][1][1][2])
        principal_axis = [round(p / max(principal_axis, key=abs), 3) for p in principal_axis]

    plt.text(0,cb+1.30, "Point group: "+fancy_subscript(Pointgroup, False), fontsize=16)
    plt.text(0,cb+0.95, "Principal axis: "+str(principal_axis), fontsize=16)

    labels = []
    labels = plot_transitions(allowed_tr1, 1, labels)
    labels = plot_transitions(allowed_tr2, 2, labels)
    handles = generate_handles(labels)

    plot_eigen(band_info1, 1)
    plot_eigen(band_info2, 2)

    plt.xlim([-1, 12])
    plt.ylim([vb-1, cb+1])
    plt.ylabel("Eigenvalue [eV]", fontsize=15)

    #ytic_pos = [vb+i for i in range(7)]
    #ytics = list(range(7))
    #plt.yticks(ytic_pos, ytics)

    plt.xticks([])
    plt.box(on=None)
    plt.subplots_adjust(left=0.160)
    #plt.legend(handles = handles, fontsize=15)
    plt.legend(handles = handles, bbox_to_anchor=(0.35, 0.80), fontsize=15)
    plt.savefig(folder_path+"/Tr_"+plotname+".png")
    #plt.savefig(folder_path+"/Tr_"+plotname+".svg", format='svg')


    return 0

if __name__ == "__main__":


    #plot_levels(sys.argv[1], sys.argv[2])
    #plot_levels_one_spin(sys.argv[1], sys.argv[3], sys.argv[2], "EIGENVAL")
    plot_levels_and_ipr(sys.argv[1], sys.argv[2], sys.argv[3])
