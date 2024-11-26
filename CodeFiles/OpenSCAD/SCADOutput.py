# -*- coding: utf-8 -*-


def fullOutput(fileName,b,xcBV,nSCAD,t_sk,t_sp,t_TE,alpha0L,WingShape,WingStruct,AFcoords):
    """
    Write to SCAD File:
        Scalars:
            Span
            BV location
            Skin Thickness
            alpha0L
        
        Matrices:    
            WingShape
            WingStruct
            Airfoil Coords
            US Fillet Coords
            LS Fillet Coords
        
    @author: justi
    """    
    import numpy as np
    from numpy import pi
    import os
    import shutil
    from distutils.dir_util import copy_tree
    from CodeFiles.OpenSCAD.GWingSCADScript import switchSettings
    from CodeFiles.OpenSCAD.GWingSCADScript import baseScript

    SCADPath = os.path.join('Wings', fileName,'OpenSCAD')
    try:
        os.mkdir(SCADPath)
        BOSL2oldLoc = os.path.join('CodeFiles','OpenSCAD','BOSL2')
        BOSL2newLoc = os.path.join('Wings', fileName,'OpenSCAD','BOSL2')
        copy_tree(BOSL2oldLoc, BOSL2newLoc)
    except:
        print("OpenSCAD FOLDER ALREADY EXISTS")
    
    # Open Output 
    scadFile = fileName + "Full.scad"
    newSCADfile = os.path.join(SCADPath,scadFile)
    SCADoutputFile = open(newSCADfile, "w")
    SCADoutputFile.write('b    = '+ str(b) +'; // [mm] Wing Span\n')
    SCADoutputFile.write('xcBV = '+ str(xcBV) +'; // [c] chordwise BV location\n')
    SCADoutputFile.write('t_sk = '+ str(t_sk) +'; // [mm] Skin Thickness\n')
    SCADoutputFile.write('t_sp = '+ str(t_sp) +'; // [mm] Spar Thickness\n')
    SCADoutputFile.write('t_TE = '+ str(t_TE) +'; // [mm] Trailing Edge Thickness\n')
    SCADoutputFile.write(switchSettings)
    
    
    n = len(WingShape)
    # Reduce discretization for faster OpenSCAD rendering
    dTheta = pi/(nSCAD - 1)
    yLLtip = WingShape[0,0]
    yInterp = np.zeros((nSCAD,1))
    WingShapeSCAD = np.zeros((nSCAD,3))
    for i in range(nSCAD):
        yInterp[i,0] = yLLtip*np.cos(pi - i*dTheta)
        WingShapeSCAD[i,0] = 0.5*np.cos(pi - i*dTheta)
    # yInterp[nSCAD-1,0] = 0
    # WingShapeSCAD[nSCAD-1,0] = 0
    # yInterp[:,0] = np.sort(yInterp[:,0],axis = 0)
    # WingShapeSCAD[:,0] = np.sort(WingShapeSCAD[:,0],axis = 0)
    WingShapeSCAD[:,1] = np.interp(yInterp,WingShape[:,0],WingShape[:,1])[:,0]
    WingShapeSCAD[:,2] = np.interp(yInterp,WingShape[:,0],WingShape[:,2])[:,0]
    
    nSpars = len(WingStruct[0])
    WingStructSCAD = np.zeros((nSCAD,nSpars))
    # breakpoint()
    for i in range(nSpars):
        WingStructSCAD[:,i] = np.interp(yInterp,WingShape[:,0],WingStruct[:,i])[:,0]

    
    
    # Write Wing Shape
    SCADoutputFile.write('WingShape = [')
    for j in range(nSCAD):
        np.savetxt(SCADoutputFile, [WingShapeSCAD[j][0]], fmt='%.5f', newline='', header='[', footer=',', comments='', encoding=None)
        np.savetxt(SCADoutputFile, [WingShapeSCAD[j][1]], fmt='%.5f', newline='', header='', footer=',', comments='', encoding=None)
        if j != nSCAD - 1:
            foot = '],\n'
        else:
            foot = ']];\n'
        np.savetxt(SCADoutputFile, [WingShapeSCAD[j][2]], fmt='%.5f', newline='', header='', footer=foot, comments='', encoding=None)
    
    # Write Wing Structure
    SCADoutputFile.write('WingStruct = [')
    for j in range(nSCAD):
        np.savetxt(SCADoutputFile, [WingStructSCAD[j][0]], fmt='%.5f', newline='', header='[', footer=',', comments='', encoding=None)
        for i in range(1,nSpars-1):
            np.savetxt(SCADoutputFile, [WingStructSCAD[j][i]], fmt='%.5f', newline='', header='', footer=',', comments='', encoding=None)
        if j != nSCAD - 1:
            foot = '],\n'
        else:
            foot = ']];\n'
        np.savetxt(SCADoutputFile, [WingStructSCAD[j][nSpars-1]], fmt='%.5f', newline='', header='', footer=foot, comments='', encoding=None)
    
    # Write Raw Airfoil Coordinates
    
    SCADoutputFile.write('alpha0L = ' + str(alpha0L) + '; // [rad]\n')
    nAF = len(AFcoords)
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
    
    from CodeFiles.OpenSCAD.GWingSCADScript import baseScript
    SCADoutputFile.write(baseScript)
    
    SCADoutputFile.close()
    
def shortOutput(fileName,b,xcBV,nSCAD,t_sk,t_sp,t_TE,alpha0L,WingShape,WingStruct,AFcoords):
    """
    Write to SCAD File:
        Scalars:
            Span
            BV location
            Skin Thickness
            alpha0L
        
        Matrices:    
            WingShape
            WingStruct
            Airfoil Coords
            US Fillet Coords
            LS Fillet Coords
        
    @author: justi
    """    
    import numpy as np
    from numpy import pi
    import os
    import shutil
    from distutils.dir_util import copy_tree
    
    SCADPath = os.path.join('Wings', fileName,'OpenSCAD')
    try:
        os.mkdir(SCADPath)
    except:
        print("OpenSCAD FOLDER ALREADY EXISTS")
    
    # Open Output File
    scadFile = fileName + "Short.scad"
    newSCADfile = os.path.join(SCADPath,scadFile)
    
    SCADoutputFile = open(newSCADfile, "w")
    SCADoutputFile.write('b    = '+ str(b) +'; // [mm] Wing Span\n')
    SCADoutputFile.write('BV_AF   = ['+ str(xcBV) +',1]; //wing location of BV\n')
    SCADoutputFile.write('t_sk = '+ str(t_sk) +'; // [mm] Skin Thickness\n')
    SCADoutputFile.write('t_sp = '+ str(t_sp) +'; // [mm] Spar Thickness\n')
    SCADoutputFile.write('t_TE = '+ str(t_TE) +'; // [mm] Trailing Edge Thickness\n')
    SCADoutputFile.write('wingLocation= [280,0,43]; // absolute location of BV\n')
    
    n = len(WingShape)
    # Reduce discretization for faster OpenSCAD rendering
    dTheta = pi/(nSCAD - 1)
    yLLtip = WingShape[0,0]
    yInterp = np.zeros((nSCAD,1))
    WingShapeSCAD = np.zeros((nSCAD,3))
    for i in range(nSCAD):
        yInterp[i,0] = yLLtip*np.cos(pi - i*dTheta)
        WingShapeSCAD[i,0] = 0.5*np.cos(pi - i*dTheta)
    WingShapeSCAD[:,1] = np.interp(yInterp,WingShape[:,0],WingShape[:,1])[:,0]
    WingShapeSCAD[:,2] = np.interp(yInterp,WingShape[:,0],WingShape[:,2])[:,0]
    
    # nSpars = max(len(WingStruct[:]))
    nSpars = 0
    for i in range(len(WingStruct[:])): nSpars = max(nSpars, len(WingStruct[i]))
    
    WingStructMat = np.ones((n,nSpars))
    for j in range(n):
        for i in range(len(WingStruct[j])):
            # breakpoint()
            WingStructMat[j,i] = WingStruct[j][i]
            
    WingStructSCAD = np.zeros((nSCAD,nSpars))
    print("nSpars = ", nSpars)
    for i in range(nSpars):
        # breakpoint()
        WingStructSCAD[:,i] = np.interp(yInterp,WingShape[:,0],WingStructMat[:,i])[:,0]

    
    
    # Write Wing Shape
    SCADoutputFile.write('WingShape = [')
    for j in range(nSCAD):
        np.savetxt(SCADoutputFile, [WingShapeSCAD[j][0]], fmt='%.5f', newline='', header='[', footer=',', comments='', encoding=None)
        np.savetxt(SCADoutputFile, [WingShapeSCAD[j][1]], fmt='%.5f', newline='', header='', footer=',', comments='', encoding=None)
        if j != nSCAD - 1:
            foot = '],\n'
        else:
            foot = ']];\n'
        np.savetxt(SCADoutputFile, [WingShapeSCAD[j][2]], fmt='%.5f', newline='', header='', footer=foot, comments='', encoding=None)
    
    # Write Wing Structure
    SCADoutputFile.write('WingStruct = [')
    for j in range(nSCAD):
        np.savetxt(SCADoutputFile, [WingStructSCAD[j][0]], fmt='%.5f', newline='', header='[', footer=',', comments='', encoding=None)
        for i in range(1,nSpars-1):
            np.savetxt(SCADoutputFile, [WingStructSCAD[j][i]], fmt='%.5f', newline='', header='', footer=',', comments='', encoding=None)
        if j != nSCAD - 1:
            foot = '],\n'
        else:
            foot = ']];\n'
        np.savetxt(SCADoutputFile, [WingStructSCAD[j][nSpars-1]], fmt='%.5f', newline='', header='', footer=foot, comments='', encoding=None)
    
    # Write Raw Airfoil Coordinates
    
    SCADoutputFile.write('alpha0L = ' + str(alpha0L) + '; // [rad]\n')
    nAF = len(AFcoords)
    SCADoutputFile.write('AFcoords = [')
    for i in range(nAF):
        np.savetxt(SCADoutputFile, [AFcoords[i][0]], fmt='%.5f', newline='', header='[', footer=',', comments='', encoding=None)
        if i != nAF - 1:
            foot = '],\n'
        else:
            foot = ']];\n'
        np.savetxt(SCADoutputFile, [AFcoords[i][1]], fmt='%.5f', newline='', header='', footer=foot, comments='', encoding=None)
    
    SCADoutputFile.close()