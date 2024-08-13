# -*- coding: utf-8 -*-
"""
GWing: FFF Wing Generator

@author: Justin D. Valenti


To Do:
    
    * Capture Airfoil geometry in AVL
    * Capture Torsional Moments
    * Torsional Predictions
    * Define GCode Orientation on Build Plate
    
    * [DEFERRED] Improve offset functions
    * [DEFERRED] Write code for intersecting internal structure
    * [DEFERRED] Overhang constraints
    * [DEFERRED] Interlocking tabs between Sections
    * [DEFERRED] Sacrificial Supporting Spars, print LE to TE after regular spars
    * [DEFERRED] Print internal structure laying down, separate from skin
"""

# ----------------------------------------------------------------------------#
# ---------------------------- Availible Wings -------------------------------#
# ---------------------------- (Uncomment One) -------------------------------#
# ----------------------------------------------------------------------------#

# from Wings.Type00.Settings import *
# from Wings.Type00_SparEx.Settings import *
# from Wings.Type10.Settings import *
# from Wings.Type10_cuffed.Settings import *
# from Wings.Type10Disp.Settings import *
# from Wings.Type11.Settings import *
# from Wings.Type11_assym.Settings import *
# from Wings.Type11_bell.Settings import *
# from Wings.Type11_thunder.Settings import *
# from Wings.Type20.Settings import *
# from Wings.Type20_parabolic.Settings import *
# from Wings.Type21.Settings import *
# from Wings.Type22.Settings import *

# from Wings.CLSparsEx1.Settings import *
# from Wings.CLSparsEx2.Settings import *
# from Wings.CLSparsEx2Wskin.Settings import *
# from Wings.CLSparsEx2DxConstr.Settings import *

# from Wings.SkinDef1BLwing.Settings import *

# from Wings.BendRect.Settings import *
from Wings.BendValenti.Settings import *
# from Wings.JoesWing.Settings import *
# from Wings.DrewWing.Settings import *

# from Wings.dxConstrained.Settings import *

# ----------------------------------------------------------------------------# 
# ----------------------------------------------------------------------------# 
# ----------------------------------------------------------------------------# 

# import importlib
# WingFolderPath = importlib.import_module("Wings."+WingFolder)
# from WingFolderPath.settings import *

import numpy as np
from CodeFiles.Aerodynamics.Airfoils.operations import *
import copy
import os
import time
import shutil

t0 = time.perf_counter()

print("Aspect Ratio = %8.0i" % AR)

try:
    print("Wingspan     = %8.0i mm" % b)
except NameError:
    pass

try:
    print("Wing Area    = %10.5f m^2" % S)
except NameError:
    pass




try:
    suppressPlot
    print("Plotting Suppressed")
except NameError:
    suppressPlot = 0

# if 'b' in globals():
#     S = b**2/AR
# elif 'S' in globals():
#     b = np.sqrt(AR*S)*1000
    
  


# print("Scale     = %6.4f" %scale)
# print("Wing Span = %6.4fm" %(b/1000))
# print("Wing Area = %6.4fm^2" %(S))
# print("Mass  = %6.4fkg" %(P/9.8))


# if GCodeWrite == 1 or SCADWrite != 0:
#     if  os.path.isdir(fileName) and fileOverwrite == 0:
#         folderName = fileName+"_new"
#     elif os.path.isdir(fileName) and fileOverwrite == 1:
#         import shutil
#         shutil.rmtree(fileName)
#         folderName = fileName
#     else:
#         folderName = fileName
#     os.mkdir(folderName)

W_G = copy.deepcopy(P)
try:
    SkinWarpTest
except NameError:
    SkinWarpTest = 0
    
if SkinWarpTest:
    plotAFCurve(c_skin,AFcoords,SparLoc)

# print("******************************************************************")
# print("******** AEROSTRUCTURAL ITERATION ",aeroStructIter," *********************")
# print("******************************************************************")
#%% Calculate Wing Shape
print("******************************************************************")
print("************************* AERODYNAMICS ***************************")
print("******************************************************************")
print("-- AERODYNAMIC DESIGN ROUTINE --")
if suppressPlot == 1:
    dispPlots = 0
else:
    dispPlots = 1

if WingType == 0.0:
    print("|------------------------|")
    print("| Explicit Wing Geometry |")
    print("| Given: Chord and Twist |")
    print("| Solve: Lift            |")
    print("|------------------------|\n")
    from CodeFiles.Aerodynamics.Design.SinglePoint import ChordTwist
    try:
        cL1 = cLaero; cL2 = cLstruct
    except:
        pass
    
    [WingShape,LiftDist] = ChordTwist(n,0.001*b,AR,cL1,cL2,\
                               a,rho,W_G,planform,twist,dispPlots)

elif WingType == 1.0:    
    print("|--------------------------|")
    print("| Single Point Aero Design |")
    print("| Given: Lift and Chord    |")
    print("| Solve: Twist             |")
    print("|--------------------------|\n")
    from CodeFiles.Aerodynamics.Design.SinglePoint import LiftChord
    try:
        cL1 = cLaero; cL2 = cLstruct
    except:
        pass
    [WingShape, LiftDist] = LiftChord(n,0.001*b,AR,cL1,cL2,a,rho,W_G,\
                                      liftShape,planform,dispPlots)
    
elif WingType == 1.1:    
    print("|-------------------------|")
    print("| Single Point Aero Design |")
    print("| Given: Lift and Twist    |")
    print("| Solve: Chord             |")
    print("|--------------------------|\n")
    from CodeFiles.Aerodynamics.Design.SinglePoint import LiftTwist
    [WingShape, LiftDist] = LiftTwist(n,0.001*b,AR,cL1,cL2,a,rho,W_G,\
                                      liftShape,twist,fileName,nAVL,pLL_AF[0],dispPlots)
        # n,b,AR,cL1,cL2,a,rho,W_G,liftShape,twist,fileName,nAVL,pLL_AF[0],dispPlots

elif WingType == 2.0:    
    print("|--------------------------|")
    print("| Two Point Aero Design    |")
    print("| Given: 2 Lift Dists      |")
    print("| Solve: Chord and Twist   |")
    print("|--------------------------|\n")
    from CodeFiles.Aerodynamics.Design.TwoPoint import TwoLift
    # WingShape = Lift_cl(n,0.001*b,AR,cL1,a,rho,W_G,liftShape,clShape,dispPlots)
    WingShape = TwoLift(n,0.001*b,AR,cL1,cL2,a,rho,W_G,liftShape1,liftShape2,dispPlots)
    
elif WingType == 2.1:    
    print("|---------------------------|")
    print("| Two Point Aero Design     |")
    print("| Given: 2 Lift Coeff Dists |")
    print("| Solve: Chord and Twist    |")
    print("|---------------------------|\n")
    from CodeFiles.Aerodynamics.Design.TwoPoint import TwoLiftCoefficient
    [WingShape,cL1,cL2] = TwoLiftCoefficient(n,0.001*b,AR,a,rho,W_G,clShape1,clShape2,dispPlots)

elif WingType == 2.2:    
    print("|--------------------------|")
    print("| Two Point Aero Design    |")
    print("| Given: Lift and cl       |")
    print("| Solve: Chord and Twist   |")
    print("|--------------------------|\n")
    from CodeFiles.Aerodynamics.Design.TwoPoint import Lift_cl
    try:
        cL1 = cLaero; cL2 = cLstruct
    except:
        cLaero = cL1; cLstruct = cL2
    WingShape = Lift_cl(n,0.001*b,AR,cL1,a,rho,W_G,liftShape,clShape,dispPlots)

from CodeFiles.Aerodynamics.Analysis.functions import PlotWing
PlotWing(WingShape[:,0],WingShape[:,1],WingShape[:,2],0.001*b,pLL_AF[0])
print("Wing Shape Generated\n")


print("-- AERODYNAMIC ANALYSIS --")
if AVLcompare:
    from CodeFiles.Aerodynamics.Analysis.functions import AVLcompare
    if WingType == 2.2:
        cLaero   = cL1
        # cLstruct = cL2
    # breakpoint()
    [AeroDists,AeroScalars] = AVLcompare(cLanalyze, WingShape, a, fileName, 0.001*b, nAVL, pLL_AF[0], dispPlots)
else:
    from CodeFiles.Aerodynamics.Analysis.functions import LLanalyzeCL
    cLaero = cL1
    cLstruct = cL2
    [AeroDists,AeroScalars] = LLanalyzeCL(WingShape,a,[cLaero,cLstruct,0],dispPlots)
LiftDist = AeroDists[:,0,0]

t1 = time.perf_counter()
aeroTime = t1 - t0


  
#%% Calculate Initial Forces and Moments
print("-- CALCULATE FORCES AND MOMENTS --")
from CodeFiles.Structures.functions import SkinWeightCalc
from CodeFiles.Structures.functions import Forces_n_Moments
from CodeFiles.Structures.functions import FuseWeightDistCalc

# print("Skin Weight = ",( 0.001*b/n)*sum(SkinWeightCalc(WingShape[:,1],AFcoords,0.001*b,0.001*extrWidth,rho_m,g)))
fuseWeightDist = FuseWeightDistCalc(P,0.001*b,n,fuseY,fuseW,WingShape[:,0])
if accountSkinWeight:
    skinWeightDist = SkinWeightCalc(WingShape[:,1],AFcoords,0.001*b,0.001*extrWidth,rho_m,g)
    skinWeight = ( 0.001*b/n)*sum(skinWeightDist)[0]
    print("Skin Weight = ",1000*skinWeight/g,'g')
else:
    print("Not Accounting for Skin Weight")

weightConverged     = 0
aeroStructIter      = 0
oneMoreTime         = 0

while (weightConverged == 0 and aeroStructIter <= aeroStructIterMax) or oneMoreTime:
    aeroStructIter += 1
    if suppressPlot == 0:
        if iterateWeights == 0 or oneMoreTime:
            dispPlots = 1
        else:
            dispPlots = 0
    else:
        dispPlots = 0
    if aeroStructIter == 1:
        structWeightDist = np.zeros((n,1))
    if accountSkinWeight:
        nonStructWeightDist = np.concatenate((fuseWeightDist,skinWeightDist),axis=1)
        LiftDist = (sum(fuseWeightDist + skinWeightDist + structWeightDist)/sum(LiftDist))*LiftDist
        
    else:
        nonStructWeightDist = fuseWeightDist
        LiftDist = (sum(fuseWeightDist + structWeightDist)/sum(LiftDist))*LiftDist
    
    BendingMomDist = Forces_n_Moments(0.001*b,WingShape[:,0],LiftDist,nonStructWeightDist,structWeightDist,dispPlots,tipLoaded)


    #%% Calculate Wing Structure
    print("-- STRUCTURAL DESIGN ROUTINE --")
    
    from CodeFiles.Structures.functions import RibLocations
    # breakpoint()
    sectBreaks = RibLocations(b,buildVolume[2],AileronStart)
    
    if StructureType == 0.0:
        from CodeFiles.Structures.Explicit import BasicSpars
        WingStruct = BasicSpars(WingShape,SparLoc)
    elif StructureType == 0.1:
        from CodeFiles.Structures.Explicit import BasicSpars
        SparLoc = genfromtxt(StructFile, delimiter=',') 
        WingStruct = Explicit.BasicSpars(WingShape,SparLoc)
    elif StructureType == 1.0 and SwitchAileron == 0:
        from CodeFiles.Structures.Adaptive import SimpleSpars
        [WingStruct,structWeightDist] = SimpleSpars(WingShape,BendingMomDist,
                                                    AFcoords, 0.001*b,
                                                    pLL_AF[0],sigma_y/FS,
                                                    rho_m,g,0.001*extrWidth)
    elif StructureType == 1.0 and SwitchAileron == 1:
        from CodeFiles.Structures.Adaptive import SimpleSparsWAil
        try:
            [WingStruct,structWeightDist,BendingInertia, FSroot,Ixroot,StructFig] = SimpleSparsWAil(WingShape,BendingMomDist,
                                                        AFcoords, 0.001*b,
                                                        pLL_AF[0],sigma_y,FS,
                                                        rho_m,g,0.001*extrWidth,
                                                        dxLim,dxDfAM,sectBreaks,dispPlots,skinStructConsider, UScurveOnly)
        except:
            [WingStruct,structWeightDist,BendingInertia, FSroot,Ixroot,StructFig] = SimpleSparsWAil(WingShape,BendingMomDist,
                                                        AFcoords, 0.001*b,
                                                        pLL_AF[0],sigma_y,FS,
                                                        rho_m,g,0.001*extrWidth,
                                                        dxLim,dxDfAM,sectBreaks,dispPlots,skinStructConsider)
    
    from CodeFiles.Structures.functions import WingBend    
    TipDisplacement = WingBend(0.001*b,E_flex,BendingInertia,LiftDist,nonStructWeightDist,structWeightDist,dispPlots)
    
    print("Wing Structure Generated")
    structWeight = (0.001*b/n)*sum(structWeightDist)[0]
    print("Structural Weight = ",1000*structWeight/g,'g')
    
    print("Initial Gross Weight= ", 1000*W_G/g," g")
    if accountSkinWeight:
        mTot = (P + skinWeight + structWeight)/g
        print("New Gross Weight    = ", 1000*mTot," g")
    else:
        mTot = (P + structWeight)/g
        print("New Gross Weight    = ", 1000*mTot," g")
    
    if iterateWeights == 0:
        print("Skipped Iteration! Proceeding with initial structure...")
        weightConverged = 1
    elif abs(W_G - (P + skinWeight + structWeight))/W_G < 10**-5:
        print("Lift and Weight Converged!")
        weightConverged = 1
    else:
        print("Not Converged!")
        W_G = P + skinWeight + structWeight
    
    if weightConverged == 1 and iterateWeights == 1 and oneMoreTime == 0:
        print("*--------------------------*")
        print("*--------------------------*")
        print("    ONE MORE TIME!!!!!      ")
        print("PICS OR IT DIDN'T HAPPEN!!!!")
        print("*--------------------------*")
        print("*--------------------------*")
        oneMoreTime = 1
    elif weightConverged == 1 and iterateWeights == 1 and oneMoreTime == 1:
        oneMoreTime = 0

     

if 0: # Output Scaling Results
    if FSroot[1]<=FS:
        sparsActive = 1
    else:
        sparsActive = 0
    if FSroot[0]<=FS:
        wingFail = 1
    else:
        wingFail = 0
        
    print("Scale           = %8.4f" %scale)
    print("Aspect Ratio    = %5.1f" %AR)
    print("Taper Ratio     = %5.1f" %planform)
    print("Wing Span       = %8.4fm" %(b/1000))
    print("Wing Area       = %8.4fm^2" %(S))
    print("Non-Wing Mass   = %8.4fkg" %(P/9.8))
    print("GTO Mass        = %8.4fkg" %(mTot))
    print("AF Thickness    = %2.0i" %(AFthickness*100))
    print("FS design       = %5.1f" %(FS))
    print("FS @ Root       = %5.1f" %(FSroot[0]))
    print("FS @ Root,Skin  = %5.1f" %(FSroot[1]))
    print("FS @ Root,Spars = %5.1f" %(FSroot[2]))
    print("delta_t/b @ 1g  = %5.4e mm " %(1000*TipDisplacement[0]/b))
    print("delta_t/b @ 3g  = %5.4e mm " %(1000*TipDisplacement[1]/b))
    print("delta_t/b @ 5g  = %5.4e mm " %(1000*TipDisplacement[2]/b))
    print("Spars Active    = %3.0i" % sparsActive)
    print("Wing Fail       = %3.0i" % wingFail)
    
    
    ScalingNumbers = [scale,AR,b/1000,S,P/9.8,mTot,AFthickness,FS,FSroot[0],FSroot[1],FSroot[2],sparsActive]
    scalingResults = open("ScalingStudy\\ScalingResults.csv", "a")  # append mode
    scalingResults.write("%8.4f," %scale)
    scalingResults.write("%5.1f," %AR)
    scalingResults.write("%3.2f," %planform)
    scalingResults.write("%8.4f," %(b/1000))
    scalingResults.write("%8.4f," %(S))
    scalingResults.write("%8.4f," %(P/9.8))
    scalingResults.write("%8.4f," %(mTot))
    scalingResults.write("%3.2f," %(AFthickness))
    scalingResults.write("%6.3f," %(extrWidth))
    scalingResults.write("%6.4e," %(FSroot[3]))
    scalingResults.write("%5.1f," %(FS))
    scalingResults.write("%5.1f," %(FSroot[0]))
    scalingResults.write("%5.1f," %(FSroot[1]))
    scalingResults.write("%5.1f," %(FSroot[2]))
    scalingResults.write("%8.4e," %(Ixroot[0]))
    scalingResults.write("%8.4e," %(Ixroot[1]))
    scalingResults.write("%8.4e," %(Ixroot[2]))
    scalingResults.write("%5.4e," %(1000*TipDisplacement[0]/b))
    scalingResults.write("%5.4e," %(1000*TipDisplacement[1]/b))
    scalingResults.write("%5.4e," %(1000*TipDisplacement[2]/b))
    scalingResults.write("%3.0i," % sparsActive)
    scalingResults.write("%3.0i," % wingFail)
    scalingResults.write("%3.0i," % (sparsActive + wingFail))
    scalingResults.write("%3.0i" % (skinScale))
    scalingResults.close()

    StructFig.savefig('ScalingStudy\image\StructFig.png')

t2 = time.perf_counter()
strucTime = t2 - t1 

#%% GCode Generation

if GCodeWrite == 1:
    print("******************************************************************")
    print("********************** WRITING GCODE *****************************")
    print("******************************************************************")
    
    try:
        flangeLength
    except:
        flangeLength = 0
        
    folderName = "Wings\\"+fileName+"\\AVL\\"
    
    from CodeFiles.PathPlanning.WingWriter import*
    GCodePath = os.path.join('Wings', fileName,'GCode')

    try:
        os.mkdir(GCodePath)
    except OSError as error:
        print("OVERWRITING EXISTING GCODE FILES")
        shutil.rmtree(GCodePath)
        os.mkdir(GCodePath)

    fileGCodeOut = os.path.join(GCodePath, fileName)
   
    from CodeFiles.Aerodynamics.Airfoils.operations import AFGeomAnalysis
    AFGeom = AFGeomAnalysis(AFcoords,101)
    tauDist = np.zeros((101,2))
    tauDist[:,0] = AFGeom[:,0]
    tauDist[:,1] = AFGeom[:,3]
    xtmax   = tauDist[np.argmax(tauDist[:,1]),0]
    
    # chunk wing into printable sections and iterate through sections
    # Make Center Section --------------------------------------------------- #
    if AileronStart != 0: # Center Section
        print("Center Section: ")
        # Generate custom start and end scripts
        GCodeStart = \
            WriteGCodeStart(fileGCodeOut,tempBedLyr0, tempExtSoft, tempExtLyr0)
        GCodeEnd   = WriteGCodeEnd(tempBedHld, tempExtEnd, tempBedEnd)
        # Open a new gcode file
        fo = open(fileGCodeOut+"C.gcode", "w")
        fo.write(GCodeStart)  # write custom start script
        fo.write("\n\n")
        # calc span of and number of layers in section
        bSect = 2*AileronStart
        nLyrs = int(np.ceil((bSect - layer0Height)/layerHeight)) + 1 
        y0 = AileronStart #starting y coord of section
        # Searching through span to find initial chord and twist for section
        chord = b*np.interp(y0,b*WingShape[:,0],WingShape[:,1])
        alpha = np.interp(y0,b*WingShape[:,0],WingShape[:,2])+ (pi/180)*alpha0L
        print("chord = %.1f mm" % (chord))
        print("alpha = %.1f deg" % (alpha*180/pi))
        # Setting print location on bed
        pLL_bed = np.array([0.5*buildVolume[0]- chord*(0.5 - pLL_AF[0]),\
                                 0.5*buildVolume[1]])
        
        # Setting Z and E gcode commands
        Z = layer0Height + Zoffset
        E = skirt(fo,AFcoords,chord,alpha,pLL_AF,pLL_bed,extrWidth,speedSolInfLyr0,\
                  layer0Height,Ehop,Zhop,Zoffset,skirtOffset,skirtnum,FlowMult,\
                  SldFillMult,FilaDiameter,-1)
        # print("E = ",E)
        print("Skirt Written!")
        
        for j in range(n-1):
            if y0 >= b*WingShape[j,0] and y0 < b*WingShape[j+1,0]:
                j0 = j 
                break
        # nSpars = len(WingStruct[j0+1])
        
        for iLyrs in range(nLyrs):  # iterate through layers
            y0 = y0 - layerHeight # set span coord for layer
            
            # # interpolate chord, alpha, and structure
            chord = b*np.interp(y0,b*WingShape[:,0],WingShape[:,1])
            alpha = np.interp(y0,b*WingShape[:,0],WingShape[:,2])+ (pi/180)*alpha0L
            
            # xcStruct = np.zeros((nSpars))
            # interpolate spar positions
            for j in range(n):
                if b*WingShape[j,0] <= y0 and b*WingShape[j+1,0] > y0:
                    if len(WingStruct[j]) == len(WingStruct[j+1]):
                        structj  = WingStruct[j]
                        structj1 = WingStruct[j+1]
                    elif len(WingStruct[j]) < len(WingStruct[j+1]):
                        nSpars = len(WingStruct[j])
                        iMainj = WingStruct[j].index(xtmax)
                        iMainj1 = WingStruct[j+1].index(xtmax)
                        idif = iMainj1 - iMainj
                        structj  = WingStruct[j]
                        structj1 = WingStruct[j+1][idif:nSpars + idif]
                    else:
                        nSpars = len(WingStruct[j+1])
                        iMainj = WingStruct[j].index(xtmax)
                        iMainj1 = WingStruct[j+1].index(xtmax)
                        idif = iMainj - iMainj1
                        structj  = WingStruct[j][idif:nSpars + idif]
                        structj1 = WingStruct[j+1]
                        
                    xcStruct = ((b*WingShape[j+1,0] - y0)/(b/n))*np.array(structj) + ((y0 - b*WingShape[j,0])/(b/n))*np.array(structj1)
                    break
            # for k in range(nSpars): 
            #     breakpoint()
            #     xcStruct[k] = np.interp(y0,b*WingShape[:,0],WingStruct[:,k])
            # Trim extraneous spars at LE and TE
            xcStruct = xcStruct[xcStruct >= extrWidth/chord]     # LE
            xcStruct = xcStruct[xcStruct <= 1 - extrWidth/chord] # TE
                
            if iLyrs == 1:
                cmndBed = "\nM140 S%i\n" % (tempBedLyr)
                cmndNoz = "M104 S%i\n" % (tempExt)
                fo.write(cmndBed + cmndNoz)
             
            Z = layer0Height + Zoffset + iLyrs*layerHeight
            
            if ribSwitch != 0 and iLyrs < ribSwitch: # rib layer 0
                E = Rib(fo,iLyrs,nLyrs,AFcoords,chord,alpha,\
                              pLL_AF,pLL_bed,extrWidth,speedPeriLyr0,\
                                  speedSolInfLyr0,speedTravel,layer0Height,layerHeight,\
                                      E,Z,Ehop,flangeLength,Zhop,FlowMult,FrstLyrMult,\
                                          SldFillMult,FilaDiameter,\
                                              infillOverlap,-1)
            elif ribSwitch != 0 and iLyrs == nLyrs-1 and SkinWarpTest == 0: # rib last layer
                E = Rib(fo,iLyrs,nLyrs,AFcoords,chord,alpha,\
                              pLL_AF,pLL_bed,extrWidth,speedPerimeter,\
                                  speedSolInfill,speedTravel,layer0Height,layerHeight,\
                                      E,Z,Ehop,flangeLength,Zhop,FlowMult,FrstLyrMult,\
                                          SldFillMult,FilaDiameter,\
                                              infillOverlap,-1)
            elif ribSwitch == 0 and iLyrs == 0: # no rib layer 0
                E = Section(fo,iLyrs,nLyrs,AFcoords,xcStruct,chord,alpha,\
                          pLL_AF,pLL_bed,extrWidth,speedPeri,speedInfLyr0,\
                              speedTravel,layer0Height,layerHeight,E,Z,\
                                  Ehop,flangeLength,Zhop,FlowMult,FrstLyrMult,\
                                      FilaDiameter,infillOverlap,-1)
            else: # no rib
                E = Section(fo,iLyrs,nLyrs,AFcoords,xcStruct,chord,alpha,\
                          pLL_AF,pLL_bed,extrWidth,speedPerimeter,speedInfill,\
                              speedTravel,layer0Height,layerHeight,E,Z,\
                                  Ehop,flangeLength,Zhop,FlowMult,FrstLyrMult,\
                                      FilaDiameter,infillOverlap,-1)
            if np.mod(iLyrs,25) == 0:
                dispString = \
                    "Center Section:  %.0f/%.0f Layers Written!" \
                    % (iLyrs,nLyrs)
                print(dispString)
                
        fo.write(GCodeEnd)
        fo.write("\n\n\n")
        FilaLength  = " %.3fm" % (E/1000)
        PartMassStr = " %.1fg" % \
            (FilaDensity*(pi*(0.5*FilaDiameter/1000)**2)*(E))
        PartCostStr = " $%.2f" % \
            (FilaDensity*(pi*(0.5*FilaDiameter/1000)**2)*(E/1000)*FilaCost)
        fo.write("\n;Material Used:")
        fo.write("\n\n;"+FilaLength)
        fo.write("\n;"+PartMassStr)
        fo.write("\n;"+PartCostStr)
        
        layerHeightStr  = " Layer Height    = %.2fmm" % (layerHeight)
        layer0HeightStr = " Layer0 Height   = %.2fmm" % (layer0Height)
        extrWidthStr    = " Extrusion Width = %.2fmm" % (extrWidth)
        fo.write("\n\n\n;Print Settings:")
        fo.write("\n\n;"+layerHeightStr)
        fo.write("\n;"+extrWidthStr)
        fo.write("\n;"+layer0HeightStr)
        fo.close() # Close opened file
    
    
    if SkinWarpTest == 0:
        # Make Right and Left Wing Sections ------------------------------------- #    
        nSemiSections = int(np.ceil((b/2 - AileronStart)/buildVolume[2]))
        print("nSemiSections = ",nSemiSections)
        for iWing in [1,-1]: # Print right wing then left wing
            # iWing =  1: Right Wing
            # iWing = -1: Left Wing
            
            for iSect in range(nSemiSections): # Iterate through sections
            # for iSect in range(1): # Test on 1 section
                if iWing == 1:
                    print("Right Section: ",iSect+1,"/",nSemiSections)
                else:
                    print("Left Section: ",iSect+1,"/",nSemiSections)
                # Generate custom start and end scripts
                GCodeStart = \
                    WriteGCodeStart(fileGCodeOut,tempBedLyr0, tempExtSoft, tempExtLyr0)
                GCodeEnd   = WriteGCodeEnd(tempBedHld, tempExtEnd, tempBedEnd)
                
                # Open a new gcode file
                if iWing == 1: #
                    fo = open(fileGCodeOut+"R"+str(iSect)+".gcode", "w")
                else:
                    fo = open(fileGCodeOut+"L"+str(iSect)+".gcode", "w")
                fo.write(GCodeStart)  # write custom start script
                fo.write("\n\n")
                
                # calc span of and number of layers in section
                bSect = (b/2 - AileronStart)/nSemiSections 
                nLyrs = int(np.ceil((bSect - layer0Height)/layerHeight)) + 1 
                y0 = iWing*((iSect/nSemiSections)*(b/2 - AileronStart) + AileronStart) # calc starting y coord of section
                
                # Searching through span to find initial chord and twist for section
                chord = b*np.interp(y0,b*WingShape[:,0],WingShape[:,1])
                alpha = np.interp(y0,b*WingShape[:,0],WingShape[:,2])+ (pi/180)*alpha0L
                
                print("chord = %.1f mm" % (chord))
                print("alpha = %.1f deg" % (alpha*180/pi))
                
                # Setting print location on bed
                if iWing == 1:
                    pLL_bed = np.array([0.5*buildVolume[0]+ chord*(0.5 - pLL_AF[0]),\
                                     0.5*buildVolume[1]])
                else:
                    pLL_bed = np.array([0.5*buildVolume[0]- chord*(0.5 - pLL_AF[0]),\
                                     0.5*buildVolume[1]])
                    
                # Setting Z and E gcode commands
                Z = layer0Height + Zoffset
                E = skirt(fo,AFcoords,chord,alpha,pLL_AF,pLL_bed,extrWidth,speedSolInfLyr0,\
                          layer0Height,Ehop,Zhop,Zoffset,skirtOffset,skirtnum,FlowMult,\
                          SldFillMult,FilaDiameter,iWing)
                # print("E = ",E)
                print("Skirt Written!")
                
                
                for j in range(n-1):
                    if y0 >= b*WingShape[j,0] and y0 < b*WingShape[j+1,0]:
                        j0 = j 
                        break
                if iWing == 1:
                    nSpars = len(WingStruct[j0])
                else:
                    nSpars = len(WingStruct[j0+1])
                
                for iLyrs in range(nLyrs):  # iterate through layers
                    y0 = y0 + iWing*layerHeight # set span coord for layer
                    
                    # # interpolate chord, alpha, and structure
                    chord = b*np.interp(y0,b*WingShape[:,0],WingShape[:,1])
                    alpha = np.interp(y0,b*WingShape[:,0],WingShape[:,2])+ (pi/180)*alpha0L
                    
                    for j in range(n - 1):
                        if b*WingShape[j,0] <= y0 and b*WingShape[j+1,0] > y0:
                            if len(WingStruct[j]) == len(WingStruct[j+1]):
                                structj  = WingStruct[j]
                                structj1 = WingStruct[j+1]
                            elif len(WingStruct[j]) < len(WingStruct[j+1]):
                                nSpars = len(WingStruct[j])
                                iMainj = WingStruct[j].index(xtmax)
                                iMainj1 = WingStruct[j+1].index(xtmax)
                                idif = iMainj1 - iMainj
                                structj  = WingStruct[j]
                                structj1 = WingStruct[j+1][idif:nSpars+idif]
                            else:
                                nSpars = len(WingStruct[j+1])
                                iMainj = WingStruct[j].index(xtmax)
                                iMainj1 = WingStruct[j+1].index(xtmax)
                                idif  = iMainj - iMainj1
                                structj  = WingStruct[j][idif:nSpars+idif]
                                structj1 = WingStruct[j+1]
                            
                            try:
                                xcStruct = ((b*WingShape[j+1,0] - y0)/(b/n))*np.array(structj) + ((y0 - b*WingShape[j,0])/(b/n))*np.array(structj1)
                            except:
                                breakpoint()
                            break
                    # Trim extraneous spars at LE and TE
                    xcStruct = xcStruct[xcStruct >= extrWidth/chord]     # LE
                    xcStruct = xcStruct[xcStruct <= 1 - extrWidth/chord] # TE
                        
                    if iLyrs == 1:
                        cmndBed = "\nM140 S%i\n" % (tempBedLyr)
                        cmndNoz = "M104 S%i\n" % (tempExt)
                        fo.write(cmndBed + cmndNoz)
                     
                    Z = layer0Height + Zoffset + iLyrs*layerHeight
                    # Case for Simple Rib, layer 0
                    if SwitchAileron == 0 and ribSwitch != 0 and iLyrs < ribSwitch:
                        E = Rib(fo,iLyrs,nLyrs,AFcoords,chord,alpha,\
                                    pLL_AF,pLL_bed,extrWidth,speedPeriLyr0,\
                                        speedSolInfLyr0,speedTravel,layer0Height,\
                                            E,Z,Ehop,flangeLength,Zhop,FlowMult,FrstLyrMult,\
                                                SldFillMult,FilaDiameter,\
                                                    infillOverlap,iWing)
                    # Case for Simple rib, wing tip
                    elif SwitchAileron == 0 and iSect + 1 == nSemiSections and iLyrs + 1 == nLyrs: 
                        E = Rib(fo,iLyrs,nLyrs,AFcoords,chord,alpha,\
                                    pLL_AF,pLL_bed,extrWidth,speedPerimeter,\
                                        speedSolInfill,speedTravel,layer0Height,\
                                            E,Z,Ehop,flangeLength,Zhop,FlowMult,FrstLyrMult,\
                                                SldFillMult,FilaDiameter,\
                                                    infillOverlap,iWing)
                    # Case for Simple Section, layer = 0
                    elif SwitchAileron == 0 and iLyrs == 0 and ribSwitch == 0:
                        E = Section(fo,iLyrs,nLyrs,AFcoords,xcStruct,chord,alpha,\
                                    pLL_AF,pLL_bed,extrWidth,speedPeriLyr0,speedInfLyr0,\
                                        speedTravel,layer0Height,layerHeight,E,Z,\
                                            Ehop,flangeLength,Zhop,FlowMult,FrstLyrMult,\
                                                FilaDiameter,infillOverlap,iWing)
                    # Case for Simple Section, layer != 0
                    elif SwitchAileron == 0 and iLyrs != 0:
                        E = Section(fo,iLyrs,nLyrs,AFcoords,xcStruct,chord,alpha,\
                                    pLL_AF,pLL_bed,extrWidth,speedPerimeter,speedInfill,\
                                        speedTravel,layer0Height,layerHeight,E,Z,\
                                            Ehop,flangeLength,Zhop,FlowMult,FrstLyrMult,\
                                                FilaDiameter,infillOverlap,iWing)
                    # Case for Rib w/Aileron, layer = 0
                    elif SwitchAileron == 1 and iLyrs < ribSwitch and ribSwitch != 0:
                        E = RibWAil(fo,iLyrs,nLyrs,AFcoords,chord,alpha,\
                              pLL_AF,pLL_bed,extrWidth,speedPeriLyr0,speedSolInfLyr0,\
                                  speedTravel,layer0Height,layerHeight,\
                                      E,Z,Ehop,flangeLength,Zhop,FlowMult,FrstLyrMult,SldFillMult,\
                                          FilaDiameter,infillOverlap,iWing,HingeGap,HingeOverlap)
                    # Case for Rib w/Aileron, layer != 0
                    elif SwitchAileron == 1 and iSect + 1 == nSemiSections and iLyrs + 1 == nLyrs:
                        E = RibWAil(fo,iLyrs,nLyrs,AFcoords,chord,alpha,\
                              pLL_AF,pLL_bed,extrWidth,speedPerimeter,speedSolInfill,\
                                  speedTravel,layer0Height,layerHeight,\
                                      E,Z,Ehop,flangeLength,Zhop,FlowMult,FrstLyrMult,SldFillMult,\
                                          FilaDiameter,infillOverlap,iWing,HingeGap,HingeOverlap)
                    # Case for Section w/Aileron, layer = 0
                    elif SwitchAileron == 1 and iLyrs == 0 and ribSwitch == 0:
                        E = SectionWAil(fo,iLyrs,nLyrs,AFcoords,xcStruct,chord,alpha,
                              pLL_AF,pLL_bed,extrWidth,speedPeriLyr0,speedInfLyr0,\
                                  speedTravel,layer0Height,\
                                      layerHeight,E,Z,Ehop,flangeLength,Zhop,FlowMult,FrstLyrMult,\
                                          FilaDiameter,infillOverlap,iWing,HingeGap,HingeOverlap)
                    # Case for Section w/Aileron, layer != 0
                    elif SwitchAileron == 1 and iLyrs != 0:
                        E = SectionWAil(fo,iLyrs,nLyrs,AFcoords,xcStruct,chord,alpha,
                              pLL_AF,pLL_bed,extrWidth,speedPerimeter,speedInfill,
                                  speedTravel,layer0Height,\
                                      layerHeight,E,Z,Ehop,flangeLength,Zhop,FlowMult,FrstLyrMult,\
                                          FilaDiameter,infillOverlap,iWing,HingeGap,HingeOverlap)
                    else:
                        print("LAYER ERROR")
                        breakpoint()
                    
                    if np.mod(iLyrs,50) == 0:
                        if iWing == 1:
                            dispString = \
                                "Right Section:%.0f/%.0f:  %.0f/%.0f Layers Written!" \
                                % (iSect+1,nSemiSections,iLyrs,nLyrs)
                        else:
                            dispString = \
                                "Left Section:%.0f/%.0f:  %.0f/%.0f Layers Written!" \
                                % (iSect+1,nSemiSections,iLyrs,nLyrs)
                        print(dispString)
                
                fo.write(GCodeEnd)
                fo.write("\n\n\n")
                FilaLength  = " %.3fm" % (E/1000)
                PartMassStr = " %.1fg" % \
                    (FilaDensity*(pi*(0.5*FilaDiameter/1000)**2)*(E))
                PartCostStr = " $%.2f" % \
                    (FilaDensity*(pi*(0.5*FilaDiameter/1000)**2)*(E/1000)*FilaCost)
                fo.write("\n;Material Used:")
                fo.write("\n\n;"+FilaLength)
                fo.write("\n;"+PartMassStr)
                fo.write("\n;"+PartCostStr)
                
                layerHeightStr  = " Layer Height    = %.2fmm" % (layerHeight)
                layer0HeightStr = " Layer0 Height   = %.2fmm" % (layer0Height)
                extrWidthStr    = " Extrusion Width = %.2fmm" % (extrWidth)
                fo.write("\n\n\n;Print Settings:")
                fo.write("\n\n;"+layerHeightStr)
                fo.write("\n;"+extrWidthStr)
                fo.write("\n;"+layer0HeightStr)
                
                
                
                
                fo.close() # Close opened file
    print("GCode Written!")
else:
    print("Did NOT Write GCode!")

t3 = time.perf_counter()
gcodeTime = t3 - t2

if SCADWrite == 1 or SCADWrite == 3:
    print("Full OpenSCAD file WRITTEN!")
    from CodeFiles.OpenSCAD.SCADOutput import fullOutput
    fullOutput(fileName,b,pLL_AF[0],nSCAD,extrWidth,extrWidth,2*extrWidth,alpha0L,WingShape,WingStruct,AFcoords)
if SCADWrite == 2 or SCADWrite == 3:
    print("Short OpenSCAD file WRITTEN!")
    from CodeFiles.OpenSCAD.SCADOutput import shortOutput
    shortOutput(fileName,b,pLL_AF[0],nSCAD,extrWidth,extrWidth,2*extrWidth,alpha0L,WingShape,WingStruct,AFcoords)
    
print("Done!!! :)")
t4 = time.perf_counter()
scadTime = t4 - t3
totalTime = t4 - t0

print("Run Time:")
print("Aerodynamic Design = %8.4f sec" % aeroTime)
print("Structural Design  = %8.4f sec" % strucTime)
print("GCode Generation   = %8.4f sec" % gcodeTime)
print("CAD Generation     = %8.4f sec" % scadTime)
print("Total Time         = %8.4f sec" % totalTime)

    