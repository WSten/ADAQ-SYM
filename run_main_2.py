#!/usr/bin/env python

from analysis import *
from extract import *
import sys

print("Starting symmetry analysis for a defect!\n")

wf_file = "WAVECAR"
HOB=find_HOB()

iprs1 = calc_ipr(wf_file,1,HOB)
ipr1_avg = np.average(iprs1)
s1bands = []
for n, ipr in enumerate(iprs1):
    if ipr > 2*ipr1_avg:
        print(HOB-15+n, ipr)
        s1bands.append(HOB-15+n)

iprs2 = calc_ipr(wf_file,2,HOB)
ipr2_avg = np.average(iprs2)
s2bands = []
for n, ipr in enumerate(iprs2):
    if ipr > 2*ipr2_avg:
        #print(HOB-15+n, ipr)
        s2bands.append(HOB-15+n)

main(s1bands, s2bands)

print("Done with symmetry analysis for a defect!")
