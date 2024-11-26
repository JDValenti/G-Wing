# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 12:16:56 2021

@author: justi
"""

import numpy as np
from numpy import pi

# import airfoil
import os
absPath = os.getcwd()
GWingPath = absPath[0:absPath.find("GWing")]+'GWing'
os.chdir(GWingPath)
from CodeFiles.Aerodynamics.Airfoils.AG24 import *
nAF = len(AFcoords)


# Write AF Coords to new file
SCADoutputFile = open("CodeFiles/OpenSCAD/Airfoils/AG24_OpenSCAD.scad", "w")
SCADoutputFile.write('AFcoords = [')
for i in range(nAF):
    np.savetxt(SCADoutputFile, [AFcoords[i][0]], fmt='%.5f', newline='', header='[', footer=',', comments='', encoding=None)
    if i != nAF - 1:
        foot = '],\n'
    else:
        foot = ']];\n'
    np.savetxt(SCADoutputFile, [AFcoords[i][1]], fmt='%.5f', newline='', header='', footer=foot, comments='', encoding=None)

# Analyze Airfoil
from CodeFiles.Aerodynamics.Airfoils.operations import AFGeomAnalysis
nInterp = 40
AFGeom = AFGeomAnalysis(AFcoords,nInterp)

# Write US points
SCADoutputFile.write('\nUScoords = [')
for i in range(nInterp):
    np.savetxt(SCADoutputFile, [AFGeom[i][0]], fmt='%.5f', newline='', header='[', footer=',', comments='', encoding=None)
    if i != nInterp - 1:
        foot = '],\n'
    else:
        foot = ']];\n'
    np.savetxt(SCADoutputFile, [AFGeom[i][1]], fmt='%.5f', newline='', header='', footer=foot, comments='', encoding=None)

# Write LS points
SCADoutputFile.write('\nLScoords = [')
for i in range(nInterp):
    np.savetxt(SCADoutputFile, [AFGeom[i][0]], fmt='%.5f', newline='', header='[', footer=',', comments='', encoding=None)
    if i != nInterp - 1:
        foot = '],\n'
    else:
        foot = ']];\n'
    np.savetxt(SCADoutputFile, [AFGeom[i][2]], fmt='%.5f', newline='', header='', footer=foot, comments='', encoding=None)

# generate and output US "fillet"
SCADoutputFile.write('\nUSfillet = [')
for i in range(nInterp):
    np.savetxt(SCADoutputFile, [AFGeom[nInterp - 1 - i][0]], fmt='%.5f', newline='', header='[', footer=',', comments='', encoding=None)
    foot = '],\n'
    np.savetxt(SCADoutputFile, [AFGeom[nInterp - 1 - i][1]], fmt='%.5f', newline='', header='', footer=foot, comments='', encoding=None)
for i in range(1,nInterp):
    np.savetxt(SCADoutputFile, [AFGeom[i][0]], fmt='%.5f', newline='', header='[', footer=',', comments='', encoding=None)
    if i != nInterp - 1:
        foot = '],\n'
    else:
        foot = ']];\n'
    np.savetxt(SCADoutputFile, [AFGeom[i][4]], fmt='%.5f', newline='', header='', footer=foot, comments='', encoding=None)

# generate and output LS "fillet"    
SCADoutputFile.write('\nLSfillet = [')
for i in range(nInterp):
    np.savetxt(SCADoutputFile, [AFGeom[nInterp - 1 - i][0]], fmt='%.5f', newline='', header='[', footer=',', comments='', encoding=None)
    foot = '],\n'
    np.savetxt(SCADoutputFile, [AFGeom[nInterp - 1 - i][4]], fmt='%.5f', newline='', header='', footer=foot, comments='', encoding=None)
for i in range(1,nInterp):
    np.savetxt(SCADoutputFile, [AFGeom[i][0]], fmt='%.5f', newline='', header='[', footer=',', comments='', encoding=None)
    if i != nInterp - 1:
        foot = '],\n'
    else:
        foot = ']];\n'
    np.savetxt(SCADoutputFile, [AFGeom[i][2]], fmt='%.5f', newline='', header='', footer=foot, comments='', encoding=None)

SCADoutputFile.close()