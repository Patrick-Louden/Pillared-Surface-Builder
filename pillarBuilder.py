#!/usr/bin/env python  

#Mark5 : Generates a double sided FCC lattice with independent pillar and inter-pillar spacing whereas Mark4 generated a double sided FCC lattice with equal-distant pillar and inter-pillar spacing.

import sys
import math
import numpy as np

#Below values for sigma and epsilon come from 
#PNAS vol. 106 no. 21 pg. 8435 (2009)   
sig = 3.4 #Angstroms
eps = 0.2325 #kJ/mol
a = sig*math.pow(2,1./6.)*math.pow(2,1./2.)
#print a
#n[i] denotes the number of unit cells in the i dimension
n = [12, 12, 14]
nmol = n[0]*n[1]*n[2]*4
boxl = [n[0], n[1], n[2]]
pillar = [3, 3, 5]
spacing = [3, 3, 5]

#r[i][j] denotes the j-dimensional coordinate of the i-th atom 
r = [[0 for x in xrange(3)]for y in xrange(nmol)]
v = [[0 for x in xrange(3)]for y in xrange(nmol)]

#placing the atoms at their sites on the perfect lattice and
#setting initial velocities to zero
iatom = 0
for i in range(0,n[0]):
    for j in range(0,n[1]):
        for k in range(0,n[2]):
            r[iatom][0] = 0.0 + i
            r[iatom][1] = 0.0 + j
            r[iatom][2] = 0.0 + k
            v[iatom][0] = 0.0
            v[iatom][1] = 0.0
            v[iatom][2] = 0.0
            iatom = iatom + 1
            r[iatom][0] = 0.0 + i
            r[iatom][1] = 0.5 + j
            r[iatom][2] = 0.5 + k
            v[iatom][0] = 0.0
            v[iatom][1] = 0.0
            v[iatom][2] = 0.0
            iatom = iatom + 1
            r[iatom][0] = 0.5 + i
            r[iatom][1] = 0.0 + j
            r[iatom][2] = 0.5 + k
            v[iatom][0] = 0.0
            v[iatom][1] = 0.0
            v[iatom][2] = 0.0
            iatom = iatom + 1
            r[iatom][0] = 0.5 + i
            r[iatom][1] = 0.5 + j
            r[iatom][2] = 0.0 + k
            v[iatom][0] = 0.0
            v[iatom][1] = 0.0
            v[iatom][2] = 0.0
            iatom = iatom + 1

#removing surface atoms to create pillars and add them to cut_block
natoms = 0
cut_block = []
xunit = pillar[0]+spacing[0]
yunit = pillar[1]+spacing[1]
for i in range(0,iatom):
####Top side of the block####
#if the atoms are located within the pillar height in the crystal
    if r[i][2] >= (boxl[2]-pillar[2]):
#if the atoms are within the spacing associated with the pillars
        if ((r[i][0] % xunit) < pillar[0]):
            if ((r[i][1] % yunit) < pillar[1]):
#store the atoms to keep to "line" and add the line to the cut_block
                line = str(i) + "\t" + "pv" + "\t" + str(r[i][0]*a) + "\t" + str(r[i][1]*a) + "\t" + str(r[i][2]*a) + "\t" + str(v[i][0]) + "\t" + str(v[i][1]) + "\t" + str(v[i][2])
                cut_block.append(line)
                natoms = natoms + 1
#if the atoms are associated with the spaces between pillars do not store them in the cut_block 
            else:
                continue
        else:
            continue

###Bottom side of the block###
#if the atoms are on the bottom side of the perfect crystal block
    elif r[i][2] < pillar[2]:
#if the atoms are within the spacing associaated with the pillars
        if ((r[i][0] % xunit) < pillar[0]):
            if ((r[i][1] % yunit) < pillar[0]):
#store the atoms to keep to "line" and add the line to the cut_block
                line = str(i) + "\t" + "pv" + "\t" + str(r[i][0]*a) + "\t" + str(r[i][1]*a) + "\t" + str(r[i][2]*a) + "\t" + str(v[i][0]) + "\t" + str(v[i][1]) + "\t" + str(v[i][2])
                cut_block.append(line)
                natoms=natoms+1
#if the atoms are associated with the spaces between pillars do not store them in the cut_block
            else:
                continue
        else:
            continue
#if the atoms are not located near the top or bottom of the crystal, but are inside the "bulk", they should be added to the cut_block 
    else:
        line = str(i) + "\t" + "pv" + "\t" + str(r[i][0]*a) + "\t" + str(r[i][1]*a) + "\t" + str(r[i][2]*a) + "\t" + str(v[i][0]) + "\t" + str(v[i][1]) + "\t" + str(v[i][2]) 
        cut_block.append(line)
        natoms = natoms + 1


#re-numbering the atoms and writing them out to an OpenMD file
b = open("mpillarSystem.md","w")
b.write("<OpenMD version=2>" + "\n" + "  <MetaData>" + "\n" + "\n")
b.write("#include \"lj2.md\"" + "\n" + "\n")
b.write("component{" + "\n")
b.write("  " + "type = " + "\"Dummy-LJ\"" + ";" + "\n")
b.write("  " + "nMol = " + str(natoms) + ";" + "\n")
b.write("}" + "\n" + "\n")
#b.write("restraint{" + "\n")
#b.write("  restraintType = \"object\";" + "\n")
#b.write("  objectSelection = \"select Dummy-LJ\";" + "\n")
#b.write("  displacementSpringConstant = 400.0;" + "\n")
#b.write("  twistSpringConstant = 0.0;" + "\n")
#b.write("  restrainedTwistAngle = 0.0;" + "\n")
#b.write("  swingXSpringConstant = 0.0;" + "\n")
#b.write("  restrainedSwingXAngle = 0.0;" + "\n")
#b.write("  swingYSpringConstant = 0.0;" + "\n")
#b.write("  restrainedSwingYAngle = 0.0;" + "\n")
#b.write("  print = \"false\";" + "\n")
#b.write("}" + "\n")
#b.write("\n")
#b.write("useRestraints = \"true\";" + "\n")
#b.write("Restraint_file = \"cut_block.in\";" + "\n")
#b.write("\n")
b.write("ensemble = NVT;" + "\n")
b.write("forceField = \"LJ2\";" + "\n")
b.write("cutoffRadius = 9;" + "\n")
b.write("cutoffMethod = \"SHIFTED_FORCE\";" + "\n")
b.write("switchingRadius = 9;" + "\n" + "\n")
b.write("targetTemp = 225;" + "\n")
b.write("targetPressure = 1.0;" + "\n" + "\n")
b.write("tauThermostat = 1e3;" + "\n")
b.write("tauBarostat = 1e4;" + "\n" + "\n")
b.write("dt = 2;" + "\n")
b.write("runTime = 1e3;" + "\n" + "\n")
b.write("tempSet = \"true\";" + "\n")
b.write("thermalTime = 10;" + "\n")
b.write("sampleTime = 100;" + "\n")
b.write("statusTime = 10;" + "\n" + "\n")
b.write("  </MetaData>" + "\n")
b.write("  <Snapshot>" + "\n")
b.write("    <FrameData>" + "\n")
b.write("      Time: 0" + "\n")
b.write("      Hmat: {{ " + str(boxl[0]*a) +", 0, 0 }, { 0, " + str(boxl[1]*a) + ", 0 }, { 0, 0, " + str(boxl[2]*a) + " }}" + "\n")
b.write("  Thermostat: 0.0 , 0.0" + "\n")
b.write("    </FrameData>" + "\n")
b.write("    <StuntDoubles>" + "\n")
new_index = 0
for line in cut_block:
    split_line = line.split()
    split_line[0] = new_index
    split_line[1] = "pv"
    b.write("\t")
    b.write('\t'.join(map(str, split_line)))
    b.write('\n')
    new_index = new_index + 1

b.write("    </StuntDoubles>" + "\n")
b.write("  </Snapshot>" + "\n")
b.write("</OpenMD>" + "\n")
b.close()

#write out to an lj2 file to store the cut-block as a rigid body
c = open("lj2.md","w")
c.write("#ifndef __LJ2_MD__" + "\n")
c.write("#define __LJ2_MD__" + "\n")
c.write("\n")
c.write("molecule{" + "\n")
c.write("\t" + "name = \"Dummy-LJ\";" + "\n")
c.write("\n")
c.write("\t" + "atom[0]{" + "\n")
c.write("\t" + "type = \"PL\";" + "\n")
c.write("\t" + "  position( 0, 0, 0 );" + "\n")
c.write("\t" + "}" + "\n")
c.write("}" + "\n")
c.write("\n")
c.write("#endif")
