#!/usr/bin/env python

from analysis import *
import sys

print("Starting symmetry analysis for a defect!\n")

lower_b_s1 = int(sys.argv[1])
upper_b_s1 = int(sys.argv[2])
lower_b_s2 = int(sys.argv[3])
upper_b_s2 = int(sys.argv[4])

s1bands = [i for i in range(lower_b_s1, upper_b_s1+1)]
s2bands = [i for i in range(lower_b_s2, upper_b_s2+1)]

main(s1bands, s2bands)

print("Done with symmetry analysis for a defect!")
