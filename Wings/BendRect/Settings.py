from numpy import pi
# --------------------------------------------------------------------------- #
# Wing Settings
# --------------------------------------------------------------------------- #

fileName = 'BendRect'
fileOverwrite = 0 
GCodeWrite = 0        # 1 = write Gcode to File, 0 = Do not
SCADWrite  = 0       # 1 = write Full OpenSCAD File; 
                      # 2 = write short OpenSCAD File;
                      # 3 = write both OpenSCAD Files;

iterateWeights = 0    # 1 = iterate wing weight until converged, 0 = don't
aeroStructIterMax   = 10

n       = 400       # Spanwise Discretization
# n       = 20       # Spanwise Discretization
b       = 1220          # Wingspan [mm]
AR      = 10            # Aspect Ratio
pLL_AF  = [0.8,1.0]    # pos of LL on root wing section, [x/c,z/t]
P       = 0.3*9.8         # Non-structural Weight [N]
rho     = 1.225         # Air Density [kg/m3]

# Pmass   = 1         # (kg) Total non-structural weight
fuseW   = 0.05      # width of fuse, frac of span
fuseY   = 0.0       # location of fuse center, frac of span, center = 0, LWT = -0.5, RWT = 0.5
rho_m   = 1.24*1000 # (kg/m^3) PLA density 
sigma_y = 26e6      # (Pa) yield stress PLA
E_tens  = 2539e6   # (Pa) ZX Tensile Modulus
E_flex  = 2470e6   # (Pa) ZX Flexural Modulus
FS      = 3
g       = 9.8       # (m/s^2) gravitational acceleration

SwitchAileron = 1  # 1 = make ailerons, 0 = do not
AileronStart  = 25 # [mm] distance from wing centerline where ailerons start
HingeGap      = 2.75 # [*extrWidth]
HingeOverlap  = 0.75 # [*extrWidth]


nSCAD = 51

# --------------------------------------------------------------------------- #
# Airfoil Selection
# --------------------------------------------------------------------------- #

from CodeFiles.Aerodynamics.Airfoils.AG24 import *

# --------------------------------------------------------------------------- #
# Aerodynamic Design Settings
# --------------------------------------------------------------------------- #
"""
Select a Wing Type:
        0.0 = Given Chord, Twist
        
        1.0 = Given Lift and Chord, Solve Twist
        1.1 = Given Lift and Twist, Solve Chord
        
        2.0 = Given 2 Lift Distributions, Solve Chord and Twist [NOT IMPLEMENTED YET]
        2.1 = Given 2 cl Distributions, Solve Chord and Twist [NOT IMPLEMENTED YET]
        2.2 = Given a lift distributiona and a cl distribution, Solve Chord and Twist [IN PROGRESS]
        
"""
WingType = 0.0
AVLcompare = 0.0
nAVL  = [100,10]
clmax = 1.2
cLanalyze = [0.2, 0.4, 1.1]

if WingType == 0.0:
    planform =  1   
    twist    = 0
    cLaero   = 0.4      # cL for aerodynamic performance
    cLstruct = 1.2      # cL for structural sizing
    
if WingType == 1.0:
    cLaero   = 0.4      # cL for aerodynamic performance
    cLstruct = 1.2      # cL for structural sizing
    planform =  0.6
    # planform =  'e'
    liftShape = 'e'
    
if WingType == 1.1:
    cLaero   = 0.4      # cL for aerodynamic performance
    cLstruct = 1.2      # cL for structural sizing
    twist    =  -2
    # planform =  'e'
    liftShape = 'e'
    
if WingType == 2.2:
    cL1   = 0.4      # cL for aerodynamic performance
    cL2 = 1.1      # cL for structural sizing
    # planform =  0.6
    liftShape = 'e'
    clShape = [[-0.5,0.9],[0.0,1.2],[0.5,0.9]]

# --------------------------------------------------------------------------- #
# Structural Design Settings
# --------------------------------------------------------------------------- #

# ribSwitch = 1 # 1 = ribs, 0 = no ribs
# tipLoaded = 0 # 1 = design structure for tip loading only
# accountSkinWeight = 1
# dxLim = 30   # max spar spacing mm

ribSwitch = 1 # 1 = ribs, 0 = no ribs
tipLoaded = 0 # 1 = design structure for tip loading only
accountSkinWeight = 1
skinStructConsider = 1
dxLim = 30   # max spar spacing mm
# dxDfAM = [2e-5,60]   # adaptive spar constraint, [0] = eta*, [1] = r*
dxDfAM = [3e-5,60]   # adaptive spar constraint, [0] = eta*, [1] = r*
"""
Select a Structure Type:
        0.0 = explicitly define chordwise locations of spars
        0.1 = explicitly define chordwise locations of spars from file
        1.0 = Basic Curvilinear Spars without skin loading
"""
StructureType = 1.0

if StructureType == 0.0:
    SparLoc = [0.25, 0.6]
    
elif StructureType == 0.1:
    from numpy import genfromtxt
    StructFile = 'WingFiles/PredefinedStructures/EllipticalStarburst.csv'

# elif StructureType == 1.0:  # Basic Curilinear Spars
    
# --------------------------------------------------------------------------- #
# Printer Settings
# --------------------------------------------------------------------------- #



printerSwitch = 1

if printerSwitch == 1: # Creality CR-10 Settings ##############################
    
    buildVolume = [300,300,150] # mm
    
    nozzleDiam   = 0.4  # mm
    layerHeight  = 0.20 # mm
    # layer0Height = 0.30 # mm
    extrWidth    = 0.5 # mm
    Zoffset      =-0.04 # mm
    # layerHeight  = 0.627*nozzleDiam # mm  # Matching rectangular extrusion area 
    # extrWidth    = 1.25*nozzleDiam  # mm  # to nozzle area
    # layerHeight  = 0.666*nozzleDiam # mm  # Matching Oval extrusion area 
    # extrWidth    = 1.25*nozzleDiam  # mm  # to nozzle area
    layer0Height = layerHeight - Zoffset # mm
    
    infillOverlap = 0.6 # fraction of extrusion width (0 = none, 1 = 1 extrWidth)
    
    FlowMult     = 1.45
    FrstLyrMult  = 1.2
    SldFillMult  = 1.2
    infillMult   = 1.0
    
    skirtOffset = 0.18 #mm
    skirtnum    = 5
    
    FilaDiameter = 1.75 # mm
    FilaDensity  = 1240 # kg/m^3
    FilaCost 	 = 25# $US/kg
    
    Zhop = 0.5 # [mm] vertical travel during hop
    Ehop = 7.0 # [mm] retraction distance during hop
    
    speedPerimeter    = 40 # mm/s
    speedInfill       = 40 # mm/s
    speedSolInfill    = 60 # mm/s
    speedPeriLyr0     = 30 # mm/s
    speedInfLyr0      = 30 # mm/s
    speedSolInfLyr0   = 30 # mm/s
    speedTravel       = 80 # mm/s
    
    tempBedLyr0 =  65 # degC
    tempBedLyr  =   0 # degC
    tempBedHld  =   0 # degC
    tempBedEnd  =   0  # degC
    
    tempExtSoft = 170 # degC
    tempExtLyr0 = 220 # degC
    tempExt     = 210 # degC
    tempExtEnd  = 140 # degC
    
    def WriteGCodeStart(fileGCodeOut, tempBedLyr0, tempExtSoft, tempExtLyr0):
    
     	GCodeStart = """;GCode Generated by GWing (J.D.Valenti 2020)
;                for the Creality CR-10
M117 %s
M140 S%i ; start bed heating up
M104 S%i; soften filament
G28 X0 Y0; home all axes
G0 Z10 F5000 ; lift nozzle
M190 R%i ; wait for bed to reach printing temp
M109 R%i ; wait for extruder to reach first lyr temp
G28; home all axes
G0 Z10 F5000 ; lift nozzle
;Finished Custom Start GCode\n
 	""" % (fileGCodeOut, tempBedLyr0, tempExtSoft, tempBedLyr0, tempExtLyr0 )
     	return GCodeStart
     	
    def WriteGCodeEnd(tempBedHld, tempExtEnd, tempBedEnd):
    
     	GCodeEnd ="""\n\n;Begin Custom End GCode\n
M140 S%i ; start bed cooling
M104 S%i ; set hot end end temp
G28 X0  ; home X axis
M84     ; disable motors
 	""" % (tempBedHld, tempExtEnd)
     	return GCodeEnd

elif printerSwitch == 2: # Lulzbot Mini 2 Settings ############################

    buildVolume = [160,160,150] # mm
    
    nozzleDiam   = 0.5  # mm
    Zoffset      =-0.05 # mm
    # layerHeight  = 0.25 # mm
    # layer0Height = 0.30 # mm
    # extrWidth    = 0.55 # mm
    # layerHeight  = 0.627*nozzleDiam # mm  # Matching rectangular extrusion area 
    # extrWidth    = 1.25*nozzleDiam  # mm  # to nozzle area
    layerHeight  = 0.666*nozzleDiam # mm  # Matching Oval extrusion area 
    extrWidth    = 1.25*nozzleDiam  # mm  # to nozzle area
    layer0Height = layerHeight - Zoffset # mm
    
    
    infillOverlap = 0.6 # fraction of extrusion width (0 = none, 1 = 1 extrWidth)
    
    FlowMult     = 1.23
    FrstLyrMult  = 1.2
    SldFillMult  = 1.2
    infillMult   = 1.0
    
    skirtOffset = 0.18 #mm
    skirtnum    = 5
    
    FilaDiameter = 2.85 # mm
    FilaDensity  = 1240 # kg/m^3
    FilaCost 	 = 25# $US/kg
    
    Zhop = 0.5 # [mm] vertical travel during hop
    Ehop = 0.5 # [mm] retraction distance during hop
    
    speedPerimeter    = 40 # mm/s
    speedInfill       = 40 # mm/s
    speedSolInfill    = 60 # mm/s
    speedPeriLyr0     = 30 # mm/s
    speedInfLyr0      = 30 # mm/s
    speedSolInfLyr0   = 30 # mm/s
    speedTravel       = 80 # mm/s
    
    tempBedLyr0 =  60 # degC
    tempBedLyr  =  0 # degC
    tempBedHld  =  0 # degC
    tempBedEnd  =  0  # degC
    
    tempExtSoft = 140 # degC
    tempExtLyr0 = 220 # degC
    tempExt     = 210 # degC
    tempExtEnd  = 140 # degC
    
    def WriteGCodeStart(fileGCodeOut, tempBedLyr0, tempExtSoft, tempExtLyr0):
    
     	GCodeStart = """;GCode Generated by GWing (J.D.Valenti 2020)
;                for the Creality CR-10
M117 %s
M107
;This G-Code has been generated specifically for the LulzBot Mini 2
M75 ; Start GLCD Print Timer
G26 ; clear potential 'probe fail' condition
G21 ; set units to Millimetres
M107 ; disable fans
G90 ; absolute positioning
M82 ; set extruder to absolute mode
G92 E0 ; set extruder position to 0
M140 S%i ; start bed heating up
G28 ; home all axes
G0 X0 Y187 Z156 F200 ; move away from endstops
M104 S%i; soften filament
G1 E-15 F75 ; retract filament
G28 X0 Y0 ; home X and Y
M109 R%i ; wait for extruder to reach probe temp
M204 S300 ; set probing acceleration
G29 ; start auto-leveling sequence
M425 Z ; use measured Z backlash for compensation
M425 Z F0 ; turn off measured Z backlash compensation. (if activated in the quality settings, this command will automatically be ignored)
M204 S2000 ; restore standard acceleration
G1 X5 Y15 Z10 F5000 ; move up off last probe point
G4 S1 ; pause
M400 ; wait for moves to finish
M117 Heating... ; progress indicator message on LCD
M109 R%i ; wait for extruder to reach initial printing temp
M190 S%i ; wait for bed to reach printing temp
G1 Z2 E0 F75 ; prime tiny bit of filment into the nozzle
M117 Mini 2 Printing... ; progress indicator message on LCD
;Finished Custom Start GCode\n
     	""" % (fileGCodeOut, tempBedLyr0, tempExtSoft, tempExtSoft, tempExtLyr0, tempBedLyr0)
     	return GCodeStart
     	
    def WriteGCodeEnd(tempBedHld, tempExtEnd, tempBedEnd):
    
     	GCodeEnd ="""\n\n;Begin Custom End GCode\n
M400 ; wait for moves to finish
M140 S%i ; start bed cooling
M104 S0 ; disable hotend
M107 ; disable fans
G92 E5 ; set extruder to 5mm for retract on print end
M117 Cooling please wait ; progress indicator message on LCD
G1 X5 Y5 Z183 E0 F3000 ; move to cooling position
G1 E5 ; re-prime extruder
M190 R%i ; wait for bed to cool down to removal temp
G1 X145 F1000 ; move extruder out of the way
G1 Y175 F1000 ; present finished print
M140 S%i; keep temperature or cool down
M77 ; End GLCD Print Timer
G90 ; absolute positioning
M18 X Y E ; turn off x y and e axis
 	""" % (tempBedHld, tempBedHld,tempBedEnd)
     	return GCodeEnd

elif printerSwitch == 3: # Lulzbot Mini 2 Small Nozzle Settings ###############

    buildVolume = [160,160,150] # mm
    
    nozzleDiam   = 0.4  # mm
    Zoffset      =-0.05 # mm
    layerHeight  = 0.2 # mm
    extrWidth    = 0.45 # mm
    # layerHeight  = 0.627*nozzleDiam # mm  # Matching rectangular extrusion area 
    # extrWidth    = 1.25*nozzleDiam  # mm  # to nozzle area
    # layerHeight  = 0.666*nozzleDiam # mm  # Matching Oval extrusion area 
    # extrWidth    = 1.25*nozzleDiam  # mm  # to nozzle area
    layer0Height = layerHeight - Zoffset # mm
    
    
    infillOverlap = 0.6 # fraction of extrusion width (0 = none, 1 = 1 extrWidth)
    
    FlowMult     = 1.34
    FrstLyrMult  = 1.2
    SldFillMult  = 1.2
    infillMult   = 1.0
    
    skirtOffset = 0.18 #mm
    skirtnum    = 5
    
    FilaDiameter = 2.85 # mm
    FilaDensity  = 1240 # kg/m^3
    FilaCost 	 = 25# $US/kg
    
    Zhop = 0.5 # [mm] vertical travel during hop
    Ehop = 0.5 # [mm] retraction distance during hop
    
    speedPerimeter    = 40 # mm/s
    speedInfill       = 40 # mm/s
    speedSolInfill    = 60 # mm/s
    speedPeriLyr0     = 30 # mm/s
    speedInfLyr0      = 30 # mm/s
    speedSolInfLyr0   = 30 # mm/s
    speedTravel       = 80 # mm/s
    
    tempBedLyr0 =  60 # degC
    tempBedLyr  =  0 # degC
    tempBedHld  =  0 # degC
    tempBedEnd  =  0  # degC
    
    tempExtSoft = 140 # degC
    tempExtLyr0 = 220 # degC
    tempExt     = 210 # degC
    tempExtEnd  = 140 # degC
    
    def WriteGCodeStart(fileGCodeOut, tempBedLyr0, tempExtSoft, tempExtLyr0):
    
     	GCodeStart = """;GCode Generated by GWing (J.D.Valenti 2020)
;                for the Creality CR-10
M117 %s
M107
;This G-Code has been generated specifically for the LulzBot Mini 2
M75 ; Start GLCD Print Timer
G26 ; clear potential 'probe fail' condition
G21 ; set units to Millimetres
M107 ; disable fans
G90 ; absolute positioning
M82 ; set extruder to absolute mode
G92 E0 ; set extruder position to 0
M140 S%i ; start bed heating up
G28 ; home all axes
G0 X0 Y187 Z156 F200 ; move away from endstops
M104 S%i; soften filament
G1 E-15 F75 ; retract filament
G28 X0 Y0 ; home X and Y
M109 R%i ; wait for extruder to reach probe temp
M204 S300 ; set probing acceleration
G29 ; start auto-leveling sequence
M425 Z ; use measured Z backlash for compensation
M425 Z F0 ; turn off measured Z backlash compensation. (if activated in the quality settings, this command will automatically be ignored)
M204 S2000 ; restore standard acceleration
G1 X5 Y15 Z10 F5000 ; move up off last probe point
G4 S1 ; pause
M400 ; wait for moves to finish
M117 Heating... ; progress indicator message on LCD
M109 R%i ; wait for extruder to reach initial printing temp
M190 S%i ; wait for bed to reach printing temp
G1 Z2 E0 F75 ; prime tiny bit of filment into the nozzle
M117 Mini 2 Printing... ; progress indicator message on LCD
;Finished Custom Start GCode\n
     	""" % (fileGCodeOut, tempBedLyr0, tempExtSoft, tempExtSoft, tempExtLyr0, tempBedLyr0)
     	return GCodeStart
     	
    def WriteGCodeEnd(tempBedHld, tempExtEnd, tempBedEnd):
    
     	GCodeEnd ="""\n\n;Begin Custom End GCode\n
M400 ; wait for moves to finish
M140 S%i ; start bed cooling
M104 S0 ; disable hotend
M107 ; disable fans
G92 E5 ; set extruder to 5mm for retract on print end
M117 Cooling please wait ; progress indicator message on LCD
G1 X5 Y5 Z183 E0 F3000 ; move to cooling position
G1 E5 ; re-prime extruder
M190 R%i ; wait for bed to cool down to removal temp
G1 X145 F1000 ; move extruder out of the way
G1 Y175 F1000 ; present finished print
M140 S%i; keep temperature or cool down
M77 ; End GLCD Print Timer
G90 ; absolute positioning
M18 X Y E ; turn off x y and e axis
 	""" % (tempBedHld, tempBedHld,tempBedEnd)
     	return GCodeEnd