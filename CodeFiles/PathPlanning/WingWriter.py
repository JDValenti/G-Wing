  # -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 12:10:47 2020

@author: justi
"""

def skirt(fo,AFcoords,chord,alpha,pLL_AF,pLL_bed,extrWidth,speedExtrusion,layer0Height,Ehop,Zhop,Zoffset,skirtOffset,skirtnum,FlowMult,SldFillMult,FilaDiameter,iWing):
    """
    Writes the skirt of a 3DP wing section.
    
    Parameters
    ----------
    fo : file object
    AFcoords : nx2 Array, Selig-Style AF Coords
    chord : scalar
    alpha : TYPE
        DESCRIPTION.
    pLL_AF : TYPE
        DESCRIPTION.
    pLL_bed : TYPE
        DESCRIPTION.
    extrWidth : TYPE
        DESCRIPTION.
    speedExtrusion : TYPE
        DESCRIPTION.
    layer0Height : TYPE
        DESCRIPTION.
    Ehop : TYPE
        DESCRIPTION.
    Zhop : TYPE
        DESCRIPTION.
    Zoffset : TYPE
        DESCRIPTION.
    skirtOffset : TYPE
        DESCRIPTION.
    skirtnum : TYPE
        DESCRIPTION.
    FlowMult : TYPE
        DESCRIPTION.
    SldFillMult : TYPE
        DESCRIPTION.
    FilaDiameter : TYPE
        DESCRIPTION.

    Returns
    -------
    E : TYPE
        DESCRIPTION.

    """
    import numpy as np
    from numpy import pi
    from CodeFiles.Aerodynamics.Airfoils.operations import AFGeomAnalysis
    nperim = len(AFcoords)  # Define perimeter length

    # Define initial perimeter ------------------
    perimeter0 = np.array(AFcoords)
    # TE thickness
    
    # rediscretize AF from Selig to [x,yUS,yLS,t,h]
    AFGeom = AFGeomAnalysis(perimeter0,11)
    tauLL  = np.interp(pLL_AF[0],AFGeom[:,0],AFGeom[:,3]) # tau @LL
    yLS_LL = np.interp(pLL_AF[0],AFGeom[:,0],AFGeom[:,2]) # yLS @LL
    
    # Place LL at x-origin
    perimeter0[:,0] = perimeter0[:,0] - pLL_AF[0] 
    perimeter0[:,1] = perimeter0[:,1] - yLS_LL - tauLL*pLL_AF[1]
    
    # Scale and rotate to desired chord and angle of attack
    perimeter0 = chord*perimeter0
    mRotAlpha = np.matrix([[np.cos(alpha),-np.sin(alpha)],
                        [np.sin(alpha), np.cos(alpha)]])
    perimeter0 = perimeter0@mRotAlpha
    
    if iWing == 1: # Mirror about y-axis for Right Wing
        perimeter0[:,0] = -perimeter0[:,0]
    
    # Translate to final position
    perimeter0[:,0] = perimeter0[:,0] + pLL_bed[0]
    perimeter0[:,1] = perimeter0[:,1] + pLL_bed[1]
    
    # Skirts = np.array(np.zeros((nperim,2,skirtnum)))
    Skirts = np.zeros((nperim,2,skirtnum))
    for i2 in list(range(skirtnum)):  # iterate through skirt loops
        for i in list(range(nperim)): # iterate though perimeter points
            if i2 == 0:
                if i == 0:
                    P0 = np.array([[perimeter0[nperim-1,0]],[perimeter0[nperim-1,1]]])
                    P1 = np.array([[perimeter0[      0 ,0]],[perimeter0[      0 ,1]]])
                    P2 = np.array([[perimeter0[      1 ,0]],[perimeter0[     1 ,1]]])
                elif i == nperim - 1:
                    P0 = np.array([[perimeter0[i-1,0]],[perimeter0[i-1,1]]])
                    P1 = np.array([[perimeter0[ i ,0]],[perimeter0[ i ,1]]])
                    P2 = np.array([[perimeter0[ 0 ,0]],[perimeter0[ 0 ,1]]])
                else:
                    P0 = np.array([[perimeter0[i-1,0]],[perimeter0[i-1,1]]])
                    P1 = np.array([[perimeter0[ i ,0]],[perimeter0[ i ,1]]])
                    P2 = np.array([[perimeter0[i+1,0]],[perimeter0[i+1,1]]])
            else:
                if i == 0:
                    P0 = np.array([[Skirts[nperim-1,0,i2-1]],[Skirts[nperim-1,1,i2-1]]])
                    P1 = np.array([[Skirts[     0 ,0,i2-1]],[Skirts[     0 ,1,i2-1]]])
                    P2 = np.array([[Skirts[     1 ,0,i2-1]],[Skirts[     1 ,1,i2-1]]])
                elif i == nperim - 1:
                    P0 = np.array([[Skirts[i-1,0,i2-1]],[Skirts[i-1,1,i2-1]]])
                    P1 = np.array([[Skirts[ i ,0,i2-1]],[Skirts[ i ,1,i2-1]]])
                    P2 = np.array([[Skirts[ 0 ,0,i2-1]],[Skirts[ 0 ,1,i2-1]]])
                else:
                    P0 = np.array([[Skirts[i-1,0,i2-1]],[Skirts[i-1,1,i2-1]]])
                    P1 = np.array([[Skirts[ i ,0,i2-1]],[Skirts[ i ,1,i2-1]]])
                    P2 = np.array([[Skirts[i+1,0,i2-1]],[Skirts[i+1,1,i2-1]]])
    
            # Phi = 0*pi/180
            # MrotPhi = matrix([[cos(Phi),-sin(Phi)],
            #                   [sin(Phi), cos(Phi)]])
    
            # P0 = MrotPhi@P0
            # P1 = MrotPhi@P1
            # P2 = MrotPhi@P2
    
            t0 = P1 - P0
            t0 = t0/np.power(t0[0]**2 + t0[1]**2,0.5)
            t2 = P2 - P1
            t2 = t2/np.power(t2[0]**2 + t2[1]**2,0.5)
            
            if iWing == 1:
                Mrot90 = np.matrix([[np.cos(pi/2),np.sin(-pi/2)],
                                    [np.sin(pi/2),np.cos(-pi/2)]])
            else:
                Mrot90 = np.matrix([[np.cos(-pi/2),-np.sin(-pi/2)],
                                    [np.sin(-pi/2), np.cos(-pi/2)]])
    
            n0 = Mrot90@t0
            n2 = Mrot90@t2
    
            n1 = (n0 + n2)
            n1 = n1/np.power(n1[0]**2 + n1[1]**2,0.5)
    
            theta0 = np.arctan2(t0[1,0],t0[0,0])
            MrotTheta0 = np.matrix([[np.cos(-theta0),-np.sin(-theta0)],
                                    [np.sin(-theta0), np.cos(-theta0)]])
            dt = MrotTheta0@t2
    
            theta = np.arctan2(dt[1,0],dt[0,0])
            if i2 == 0:
                eps_mag = (0.5*extrWidth + skirtOffset)/(np.cos(theta/2))
            else:
                eps_mag =     (extrWidth)/(np.cos(theta/2))
            eps = eps_mag*n1
            P1offset = P1 + eps
    
            Skirts[i,0,i2] = P1offset[0,0]
            Skirts[i,1,i2] = P1offset[1,0]
    
    # ---------------------- Skirt/Brim GCode Composition ------------------- #
    if skirtOffset == 0:
        skirt_command = "\n\nM117Brim\n"
    else:
        skirt_command = "\n\nM117Skirt\n"
    cmndG = "G0"
    E = 0
    Z = layer0Height + Zoffset
    
    cmndF = " F%.2f" % (60*speedExtrusion)
    cmndZ = " Z%.4f" % (Z)
     
    skirt_command = skirt_command + cmndG
    skirt_command = skirt_command + cmndZ
    skirt_command = skirt_command + cmndF
    skirt_command = skirt_command + "  ;Set Z"
    
    for i2 in list(range(skirtnum)):
        skirt_command = skirt_command + "\n"
        cmndG = "G0"
        cmndX = " X%.4f" % (Skirts[0,0,skirtnum - 1 - i2])
        cmndY = " Y%.4f" % (Skirts[0,1,skirtnum - 1 - i2])
        skirt_command = skirt_command + cmndG
        skirt_command = skirt_command + cmndF
        skirt_command = skirt_command + cmndX
        skirt_command = skirt_command + cmndY
        if skirtOffset == 0:
            skirt_command = skirt_command + "  ;Brim Loop Start"
        else:
            skirt_command = skirt_command + "  ;Skirt Loop Start"
    
        for i1 in list(range(1,nperim)):
            skirt_command = skirt_command + "\n"
            cmndG = "G1"
    
            nzlXo = Skirts[i1-1,0,skirtnum - 1 - i2]
            nzlYo = Skirts[i1-1,1,skirtnum - 1 - i2]
            nzlX = Skirts[i1,0,skirtnum - 1 - i2]
            nzlY = Skirts[i1,1,skirtnum - 1 - i2]
            nzlTrvl = np.sqrt((nzlX - nzlXo)**2 + (nzlY - nzlYo)**2)
            cmndX = " X%.4f" % (nzlX)
            cmndY = " Y%.4f" % (nzlY)
    
            Eprime = extrWidth*layer0Height
            dE = SldFillMult*FlowMult*Eprime*nzlTrvl/FilaDiameter**2
            E = E + dE
            cmndE = " E%.4f" % (E)
    
            skirt_command = skirt_command + cmndG
            skirt_command = skirt_command + cmndF
            skirt_command = skirt_command + cmndE
            skirt_command = skirt_command + cmndX
            skirt_command = skirt_command + cmndY
            if skirtOffset == 0:
                skirt_command = skirt_command + "  ;Brim"
            else:
                skirt_command = skirt_command + "  ;Skirt"
                
        skirt_command = skirt_command + "\n"
        cmndG = "G1"
    
        nzlXo = Skirts[nperim-1,0,skirtnum - 1 - i2]
        nzlYo = Skirts[nperim-1,1,skirtnum - 1 - i2]
        nzlX = Skirts[0,0,skirtnum - 1 - i2]
        nzlY = Skirts[0,1,skirtnum - 1 - i2]
    
        nzlTrvl = np.sqrt((nzlX - nzlXo)**2 + (nzlY - nzlYo)**2)
        Eprime = extrWidth*layer0Height
        dE = 1.5*FlowMult*Eprime*nzlTrvl/FilaDiameter**2
        E = E + dE
        cmndE = " E%.4f" % (E)
        cmndX = " X%.4f" % (Skirts[0,0,skirtnum - 1 - i2])
        cmndY = " Y%.4f" % (Skirts[0,1,skirtnum - 1 - i2])
        skirt_command = skirt_command + cmndG
        skirt_command = skirt_command + cmndF
        skirt_command = skirt_command + cmndE
        skirt_command = skirt_command + cmndX
        skirt_command = skirt_command + cmndY
        if skirtOffset == 0:
            skirt_command = skirt_command + "  ;Brim"
        else:
            skirt_command = skirt_command + "  ;Skirt"
            
    # Lift Nozzle with retraction
    E = E - Ehop
    cmndG = "\nG1"
    cmndF = " F%.2f" % (60*speedExtrusion)
    cmndE = " E%.2f" % (E)
    cmndZ = " Z%.4f" % (Z+Zhop)
    skirt_command = skirt_command + cmndG
    skirt_command = skirt_command + cmndF
    skirt_command = skirt_command + cmndE
    skirt_command = skirt_command + cmndX
    skirt_command = skirt_command + cmndY
    skirt_command = skirt_command + cmndZ
    skirt_command = skirt_command + "  ;Retraction/Zhop"
            
    fo.write(skirt_command)
    return(E)

def Section(fo,iLyrs,nLyrs,AFcoords,xcStruct,chord,alpha,pLL_AF,pLL_bed,extrWidth,speedExtrusion,speedInfExtrusion,speedTravel,layer0Height,layerHeight,E,Z,Ehop,EhopP,Zhop,FlowMult,FrstLyrMult,FilaDiameter,infillOverlap,iWing):
    # breakpoint()
    import numpy as np
    from CodeFiles.Aerodynamics.Airfoils.operations import AFGeomAnalysis
    chord = float(chord) # ensuring chord is acutally a float
    nperim = len(AFcoords)  # Define perimeter length

    # Define initial perimeter ------------------
    perimeter = np.array(AFcoords)
    # TE thickness
    
    # rediscretize AF from Selig to [x,yUS,yLS,t,h]
    AFGeom = AFGeomAnalysis(perimeter,11)
    tauLL  = np.interp(pLL_AF[0],AFGeom[:,0],AFGeom[:,3]) # tau @LL
    yLS_LL = np.interp(pLL_AF[0],AFGeom[:,0],AFGeom[:,2]) # yLS @LL
    
    # Place LL at x-origin
    perimeter[:,0] = perimeter[:,0] - pLL_AF[0] 
    perimeter[:,1] = perimeter[:,1] - yLS_LL - tauLL*pLL_AF[1]
    
    # Scale and rotate to desired chord and angle of attack
    perimeter = chord*perimeter
    mRotAlpha = np.matrix([[np.cos(alpha),-np.sin(alpha)],
                        [np.sin(alpha), np.cos(alpha)]])
    perimeter = perimeter@mRotAlpha
    
    if iWing == 1: # Mirror about y-axis for Right Wing
        perimeter[:,0] = -perimeter[:,0]
    
    # Translate to final position
    perimeter[:,0] = perimeter[:,0] + pLL_bed[0]
    perimeter[:,1] = perimeter[:,1] + pLL_bed[1]
    
    # if iWing == 1:
    #     xLE = float(max(perimeter[:,0]))
    # else:
    #     xLE = float(min(perimeter[:,0]))
    xMin = float(min(perimeter[:,0]))
    
    # Find US and LS Coords at Spar Locations -------------------------------
    AFPlaced = AFGeomAnalysis(perimeter,101)
    nSpars = len(xcStruct)
    xcStruct.sort() # Sort Structure LE -> TE
    xcStruct = np.flip(xcStruct) # Flip struture direction to TE -> LE
    if iWing == 1:
        xcStruct = -xcStruct + 1
        
    infill = np.zeros((nSpars,3))
    # for i in range(nSpars):
    #     infill[i][0] = chord*xcStruct[i] + xMin
    infill[:,0] = chord*xcStruct + xMin
    infill[:,1] = np.interp(infill[:,0],AFPlaced[:,0],AFPlaced[:,1]) - infillOverlap*extrWidth
    infill[:,2] = np.interp(infill[:,0],AFPlaced[:,0],AFPlaced[:,2]) + infillOverlap*extrWidth
    
    # breakpoint()
    
    # -------------------- Perimeter GCode Composition ---------------------- #
    # Initialize perimeter gcode                
    perimeterR = np.flipud(perimeter) # reverse perimeter to out LS then US
    perimeter_command = "\n\nM117 L%.0f/%.0f;Z%.1f"% (iLyrs,nLyrs,Z)
    
    # Move to Perimeter start
    cmndG = "\nG0"
    cmndF = " F%.2f" % (60*speedTravel)
    cmndX = " X%.4f" % (perimeterR[0,0])
    cmndY = " Y%.4f" % (perimeterR[0,1])
    cmndZ = " Z%.4f" % (Z+Zhop)
    perimeter_command = perimeter_command + cmndG
    perimeter_command = perimeter_command + cmndF
    perimeter_command = perimeter_command + cmndX
    perimeter_command = perimeter_command + cmndY
    perimeter_command = perimeter_command + cmndZ
    perimeter_command = perimeter_command + "  ;Move to Perimeter Start"
    
    # Set Z Height
    cmndG = "\nG0"
    cmndF = " F%.2f" % (60*speedTravel)
    cmndX = " X%.4f" % (perimeterR[0,0])
    cmndY = " Y%.4f" % (perimeterR[0,1])
    cmndZ = " Z%.4f" % (Z)
    perimeter_command = perimeter_command + cmndG
    perimeter_command = perimeter_command + cmndF
    perimeter_command = perimeter_command + cmndX
    perimeter_command = perimeter_command + cmndY
    perimeter_command = perimeter_command + cmndZ
    perimeter_command = perimeter_command + "  ;Set Z Height"
    
    # Prime Nozzle
    E = E + Ehop + EhopP
    cmndG = "\nG1"
    cmndF = " F%.2f" % (60*speedExtrusion)
    cmndE = " E%.2f" % (E)
    cmndZ = " Z%.4f" % (Z)
    perimeter_command = perimeter_command + cmndG
    perimeter_command = perimeter_command + cmndF
    perimeter_command = perimeter_command + cmndE
    perimeter_command = perimeter_command + cmndX
    perimeter_command = perimeter_command + cmndY
    perimeter_command = perimeter_command + cmndZ
    perimeter_command = perimeter_command + "  ;Prime Nozzle"

    for i1 in range(1,nperim):
        cmndG = "\nG1"
        
        nzlXo = perimeterR[i1-1,0]
        nzlYo = perimeterR[i1-1,1]
        nzlX  = perimeterR[i1  ,0]
        nzlY  = perimeterR[i1  ,1]
        nzlTrvl = np.sqrt((nzlX - nzlXo)**2 + (nzlY - nzlYo)**2)
        cmndX = " X%.4f" % (nzlX)
        cmndY = " Y%.4f" % (nzlY)
        
        if iLyrs == 0:
            Eprime = extrWidth*layer0Height
            dE = FrstLyrMult*FlowMult*Eprime*nzlTrvl/FilaDiameter**2
        else:
            Eprime = extrWidth*layerHeight
            dE = FlowMult*Eprime*nzlTrvl/FilaDiameter**2
            
        E = E + dE
        cmndE = " E%.4f" % (E)
        
        perimeter_command = perimeter_command + cmndG
        perimeter_command = perimeter_command + cmndF
        perimeter_command = perimeter_command + cmndE
        perimeter_command = perimeter_command + cmndX
        perimeter_command = perimeter_command + cmndY
        perimeter_command = perimeter_command + "  ;Perimeter"
    
    cmndG = "\nG1"
    
    nzlXo = perimeterR[-1,0]
    nzlYo = perimeterR[-1,1]
    nzlX  = perimeterR[ 0,0]
    nzlY  = perimeterR[ 0,1]
    
    nzlTrvl = np.sqrt((nzlX - nzlXo)**2 + (nzlY - nzlYo)**2)
    if iLyrs == 0:
        Eprime = extrWidth*layer0Height
        dE = FrstLyrMult*FlowMult*Eprime*nzlTrvl/FilaDiameter**2
    else:
        Eprime = extrWidth*layerHeight
        dE = FlowMult*Eprime*nzlTrvl/FilaDiameter**2
        
    E = E + dE
    cmndE = " E%.4f" % (E)
    
    cmndX = " X%.4f" % (perimeterR[0,0])
    cmndY = " Y%.4f" % (perimeterR[0,1])
    perimeter_command = perimeter_command + cmndG
    perimeter_command = perimeter_command + cmndF
    perimeter_command = perimeter_command + cmndE
    perimeter_command = perimeter_command + cmndX
    perimeter_command = perimeter_command + cmndY
    perimeter_command = perimeter_command + cmndZ
    perimeter_command = perimeter_command + "  ;Perimeter End"
            
    fo.write(perimeter_command)
    
    # ---------------- Internal Structure GCode Composition ----------------- #
    infill_command = "\n"
    
    # Lift Nozzle with retraction
    E = E - Ehop
    cmndG = "\nG1"
    cmndF = " F%.2f" % (60*speedExtrusion)
    cmndE = " E%.2f" % (E)
    cmndZ = " Z%.4f" % (Z+Zhop)
    infill_command = infill_command + cmndG
    infill_command = infill_command + cmndF
    infill_command = infill_command + cmndE
    infill_command = infill_command + cmndX
    infill_command = infill_command + cmndY
    infill_command = infill_command + cmndZ
    infill_command = infill_command + "  ;Retraction/Zhop"
    
    for i in range(nSpars):
        # Travel to Spar Location
        cmndG = "\nG0"
        cmndF = " F%.2f" % (60*speedTravel)
        cmndX = " X%.4f" % (infill[i][0])
        if np.mod(i,2):
            cmndY = " Y%.4f" % (infill[i][1])
        else:
            cmndY = " Y%.4f" % (infill[i][2])
        infill_command = infill_command + cmndG
        infill_command = infill_command + cmndF
        infill_command = infill_command + cmndX
        infill_command = infill_command + cmndY
        infill_command = infill_command + cmndZ
        infill_command = infill_command + "  ;Travel to Spar"
        
        # Set Z Height
        cmndZ = " Z%.4f" % (Z)
        infill_command = infill_command + cmndG
        infill_command = infill_command + cmndF
        infill_command = infill_command + cmndX
        infill_command = infill_command + cmndY
        infill_command = infill_command + cmndZ
        infill_command = infill_command + "  ;Set Z Height"
        
        # Prime Nozzle
        E = E + Ehop + EhopP
        cmndG = "\nG1"
        cmndF = " F%.2f" % (60*speedExtrusion)
        cmndE = " E%.2f" % (E)
        infill_command = infill_command + cmndG
        infill_command = infill_command + cmndF
        infill_command = infill_command + cmndE
        infill_command = infill_command + cmndX
        infill_command = infill_command + cmndY
        infill_command = infill_command + cmndZ
        infill_command = infill_command + "  ;Prime Nozzle"
        
        # Spar Command
        nzlTrvl = np.sqrt((infill[i][1] - infill[i][2])**2)
        if iLyrs == 0:
            Eprime = extrWidth*layer0Height
            dE = FrstLyrMult*FlowMult*Eprime*nzlTrvl/FilaDiameter**2
        else:
            Eprime = extrWidth*layerHeight
            dE = FlowMult*Eprime*nzlTrvl/FilaDiameter**2
        E = E + dE
        cmndF = " F%.2f" % (60*speedInfExtrusion)
        cmndE = " E%.2f" % (E)
        cmndX = " X%.4f" % (infill[i][0])
        if np.mod(i,2):
            cmndY = " Y%.4f" % (infill[i][2])
        else:
            cmndY = " Y%.4f" % (infill[i][1])
        infill_command = infill_command + cmndG
        infill_command = infill_command + cmndF
        infill_command = infill_command + cmndE
        infill_command = infill_command + cmndX
        infill_command = infill_command + cmndY
        infill_command = infill_command + cmndZ
        infill_command = infill_command + "  ;Spar"
        
        # Lift Nozzle with retraction
        E = E - Ehop
        cmndG = "\nG1"
        cmndE = " E%.2f" % (E)
        cmndZ = " Z%.4f" % (Z+Zhop)
        infill_command = infill_command + cmndG
        infill_command = infill_command + cmndF
        infill_command = infill_command + cmndE
        infill_command = infill_command + cmndX
        infill_command = infill_command + cmndY
        infill_command = infill_command + cmndZ
        infill_command = infill_command + "  ;Retraction/Zhop"
    fo.write(infill_command)
    return E

def Rib(fo,iLyrs,nLyrs,AFcoords,chord,alpha,pLL_AF,pLL_bed,extrWidth,speedExtrusion,speedInfExtrusion,speedTravel,layer0Height,layerHeight,E,Z,Ehop,EhopP,Zhop,FlowMult,FrstLyrMult,SldFillMult,FilaDiameter,infillOverlap,iWing):
    # breakpoint()
    import numpy as np
    from CodeFiles.Aerodynamics.Airfoils.operations import AFGeomAnalysis
    chord = float(chord) # ensuring chord is acutally a float
    nperim = len(AFcoords)  # Define perimeter length

    # Define initial perimeter ------------------
    perimeter = np.array(AFcoords)
    # TE thickness
    
    # rediscretize AF from Selig to [x,yUS,yLS,t,h]
    AFGeom = AFGeomAnalysis(perimeter,11)
    tauLL  = np.interp(pLL_AF[0],AFGeom[:,0],AFGeom[:,3]) # tau @LL
    yLS_LL = np.interp(pLL_AF[0],AFGeom[:,0],AFGeom[:,2]) # yLS @LL
    
    # Place LL at x-origin
    perimeter[:,0] = perimeter[:,0] - pLL_AF[0] 
    perimeter[:,1] = perimeter[:,1] - yLS_LL - tauLL*pLL_AF[1]
    
    # Scale and rotate to desired chord and angle of attack
    perimeter = chord*perimeter
    mRotAlpha = np.matrix([[np.cos(alpha),-np.sin(alpha)],
                        [np.sin(alpha), np.cos(alpha)]])
    perimeter = perimeter@mRotAlpha
    
    if iWing == 1: # Mirror about y-axis for Right Wing
        perimeter[:,0] = -perimeter[:,0]
    
    # Translate to final position
    perimeter[:,0] = perimeter[:,0] + pLL_bed[0]
    perimeter[:,1] = perimeter[:,1] + pLL_bed[1]
    
    # if iWing == 1:
    #     xLE = float(max(perimeter[:,0]))
    # else:
    #     xLE = float(min(perimeter[:,0]))
    xMin = float(min(perimeter[:,0]))
    xMax = float(max(perimeter[:,0]))
    
    # Find US and LS Coords at Spar Locations -------------------------------
    AFPlaced = AFGeomAnalysis(perimeter,101)
    
    solidFillLEBuffer = 5
    solidFillTEBuffer = 30
    if iWing == 1:
        xLE = xMax
        xTE = xMin
    else:
        xLE = xMin
        xTE = xMax
        
    infill = np.zeros((1,3))
    infillStart = xTE + iWing*extrWidth*solidFillTEBuffer
    infill[0,0] = infillStart
    infill[0,1] = np.interp(infill[:,0],AFPlaced[:,0],AFPlaced[:,1]) - infillOverlap*extrWidth
    infill[0,2] = np.interp(infill[:,0],AFPlaced[:,0],AFPlaced[:,2]) + infillOverlap*extrWidth
    i = 0
    while 1:
        i += 1
        infill = np.append(infill,[[0,0,0]],axis = 0)
        infill[i,0] = infill[i-1,0] + iWing*extrWidth
        infill[i,1] = np.interp(infill[i,0],AFPlaced[:,0],AFPlaced[:,1]) - infillOverlap*extrWidth
        infill[i,2] = np.interp(infill[i,0],AFPlaced[:,0],AFPlaced[:,2]) + infillOverlap*extrWidth
        if iWing*(xLE - iWing*solidFillLEBuffer*extrWidth - infill[i,0]) < 0:
            break
            
    nFill = len(infill[:,0])
    # breakpoint()
    
    # -------------------- Perimeter GCode Composition ---------------------- #
    perimeterR = np.flipud(perimeter) # reverse perimeter to out LS then US
    # Initialize perimeter gcode                
    perimeter_command = "\n\nM117 L%.0f/%.0f;Z%.1f"% (iLyrs,nLyrs,Z)
    
    # Move to Perimeter start
    cmndG = "\nG0"
    cmndF = " F%.2f" % (60*speedTravel)
    cmndX = " X%.4f" % (perimeterR[0,0])
    cmndY = " Y%.4f" % (perimeterR[0,1])
    cmndZ = " Z%.4f" % (Z+Zhop)
    perimeter_command = perimeter_command + cmndG
    perimeter_command = perimeter_command + cmndF
    perimeter_command = perimeter_command + cmndX
    perimeter_command = perimeter_command + cmndY
    perimeter_command = perimeter_command + cmndZ
    perimeter_command = perimeter_command + "  ;Move to Perimeter Start"
    
    # Set Z Height
    cmndG = "\nG0"
    cmndF = " F%.2f" % (60*speedTravel)
    cmndX = " X%.4f" % (perimeterR[0,0])
    cmndY = " Y%.4f" % (perimeterR[0,1])
    cmndZ = " Z%.4f" % (Z)
    perimeter_command = perimeter_command + cmndG
    perimeter_command = perimeter_command + cmndF
    perimeter_command = perimeter_command + cmndX
    perimeter_command = perimeter_command + cmndY
    perimeter_command = perimeter_command + cmndZ
    perimeter_command = perimeter_command + "  ;Set Z Height"
    
    # Prime Nozzle
    E = E + Ehop + EhopP
    cmndG = "\nG1"
    cmndF = " F%.2f" % (60*speedExtrusion)
    cmndE = " E%.2f" % (E)
    cmndZ = " Z%.4f" % (Z)
    perimeter_command = perimeter_command + cmndG
    perimeter_command = perimeter_command + cmndF
    perimeter_command = perimeter_command + cmndE
    perimeter_command = perimeter_command + cmndX
    perimeter_command = perimeter_command + cmndY
    perimeter_command = perimeter_command + cmndZ
    perimeter_command = perimeter_command + "  ;Prime Nozzle"

    for i1 in range(1,nperim):
        cmndG = "\nG1"
        
        nzlXo = perimeterR[i1-1,0]
        nzlYo = perimeterR[i1-1,1]
        nzlX  = perimeterR[i1  ,0]
        nzlY  = perimeterR[i1  ,1]
        nzlTrvl = np.sqrt((nzlX - nzlXo)**2 + (nzlY - nzlYo)**2)
        cmndX = " X%.4f" % (nzlX)
        cmndY = " Y%.4f" % (nzlY)
        
        if iLyrs == 0:
            Eprime = extrWidth*layer0Height
            dE = FrstLyrMult*FlowMult*Eprime*nzlTrvl/FilaDiameter**2
        else:
            Eprime = extrWidth*layerHeight
            dE = FlowMult*Eprime*nzlTrvl/FilaDiameter**2
            
        E = E + dE
        cmndE = " E%.4f" % (E)
        
        perimeter_command = perimeter_command + cmndG
        perimeter_command = perimeter_command + cmndF
        perimeter_command = perimeter_command + cmndE
        perimeter_command = perimeter_command + cmndX
        perimeter_command = perimeter_command + cmndY
        perimeter_command = perimeter_command + "  ;Perimeter"
    
    cmndG = "\nG1"
    
    nzlXo = perimeterR[-1,0]
    nzlYo = perimeterR[-1,1]
    nzlX  = perimeterR[ 0,0]
    nzlY  = perimeterR[ 0,1]
    
    nzlTrvl = np.sqrt((nzlX - nzlXo)**2 + (nzlY - nzlYo)**2)
    Eprime = extrWidth*layer0Height
    dE = FrstLyrMult*FlowMult*Eprime*nzlTrvl/FilaDiameter**2
        
    E = E + dE
    cmndE = " E%.4f" % (E)
    
    cmndX = " X%.4f" % (perimeterR[0,0])
    cmndY = " Y%.4f" % (perimeterR[0,1])
    perimeter_command = perimeter_command + cmndG
    perimeter_command = perimeter_command + cmndF
    perimeter_command = perimeter_command + cmndE
    perimeter_command = perimeter_command + cmndX
    perimeter_command = perimeter_command + cmndY
    perimeter_command = perimeter_command + cmndZ
    perimeter_command = perimeter_command + "  ;Perimeter End"
            
    fo.write(perimeter_command)
    
    # ---------------- Rib (Solid Fill) GCode Composition ----------------- #
    infill_command = "\n"
    
    for i in range(nFill):
        # Travel to Spar Location
        cmndG = "\nG0"
        cmndF = " F%.2f" % (60*speedTravel)
        cmndX = " X%.4f" % (infill[i][0])
        if np.mod(i,2):
            cmndY = " Y%.4f" % (infill[i][1])
        else:
            cmndY = " Y%.4f" % (infill[i][2])
        infill_command = infill_command + cmndG
        infill_command = infill_command + cmndF
        infill_command = infill_command + cmndX
        infill_command = infill_command + cmndY
        infill_command = infill_command + cmndZ
        infill_command = infill_command + "  ;Travel to Fill"
        
        # Spar Command
        nzlTrvl = np.sqrt((infill[i][1] - infill[i][2])**2)
        Eprime = extrWidth*layer0Height
        dE = FrstLyrMult*FlowMult*Eprime*nzlTrvl/FilaDiameter**2
        E = E + dE
        cmndF = " F%.2f" % (60*speedInfExtrusion)
        cmndE = " E%.2f" % (E)
        cmndX = " X%.4f" % (infill[i][0])
        if np.mod(i,2):
            cmndY = " Y%.4f" % (infill[i][2])
        else:
            cmndY = " Y%.4f" % (infill[i][1])
        infill_command = infill_command + cmndG
        infill_command = infill_command + cmndF
        infill_command = infill_command + cmndE
        infill_command = infill_command + cmndX
        infill_command = infill_command + cmndY
        infill_command = infill_command + cmndZ
        infill_command = infill_command + "  ;Rib"
        
    # Lift Nozzle with retraction
    E = E - Ehop
    cmndG = "\nG1"
    cmndE = " E%.2f" % (E)
    cmndZ = " Z%.4f" % (Z+Zhop)
    infill_command = infill_command + cmndG
    infill_command = infill_command + cmndF
    infill_command = infill_command + cmndE
    infill_command = infill_command + cmndX
    infill_command = infill_command + cmndY
    infill_command = infill_command + cmndZ
    infill_command = infill_command + "  ;Retraction/Zhop"
    fo.write(infill_command)
    return E

def SectionWAil(fo,iLyrs,nLyrs,AFcoords,xcStruct,chord,alpha,pLL_AF,pLL_bed,extrWidth,speedExtrusion,speedInfExtrusion,speedTravel,layer0Height,layerHeight,E,Z,Ehop,flangeLength,Zhop,FlowMult,FrstLyrMult,FilaDiameter,infillOverlap,iWing,HingeGap,HingeOverlap):
    # breakpoint()
    import numpy as np
    from CodeFiles.Aerodynamics.Airfoils.operations import AFGeomAnalysis
    chord = float(chord) # ensuring chord is acutally a float
    nperim = len(AFcoords)  # Define perimeter length
    
    
    # Define initial perimeter ------------------
    perimeter = np.array(AFcoords)
    
    # TE thickness
    
    # rediscretize AF from Selig to [x,yUS,yLS,t,h]
    AFGeom = AFGeomAnalysis(perimeter,201)
    tauLL  = np.interp(pLL_AF[0],AFGeom[:,0],AFGeom[:,3]) # tau @LL
    yLS_LL = np.interp(pLL_AF[0],AFGeom[:,0],AFGeom[:,2]) # yLS @LL
    
    # generate airfoil with hinge
    perimeterWHinge = np.zeros((2*nperim,2))
    jHingePoints = np.zeros(4)
    HingeGap = HingeGap*extrWidth/chord
    cutWidth   = HingeGap + tauLL
    cutMade = 0
    CheckUS = 1
    CheckLS = 0
    i = -1 # initialize counter for base perimeter
    j = -1 # initialize counter for perimeter with hinge
    while i < nperim - 1:
        i += 1
        if i >2:
            CheckUS = perimeter[i,0] <  perimeter[i-1,0]
            CheckLS = perimeter[i,0] >= perimeter[i-1,0]
        
        CheckUScase1 = perimeter[i,0]   >= pLL_AF[0] + HingeGap
        CheckUScase2 = perimeter[i,0]   <  pLL_AF[0] \
                   and perimeter[i-1,0] >= pLL_AF[0] + HingeGap
        CheckUScase3 = perimeter[i,0]   >= pLL_AF[0] \
                   and perimeter[i,0]   <  pLL_AF[0] + HingeGap \
                   and perimeter[i-1,0] >= pLL_AF[0] + HingeGap
        CheckUScase3 = perimeter[i,0]   >= pLL_AF[0] \
                   and perimeter[i,0]   <  pLL_AF[0] + HingeGap \
                   and perimeter[i-1,0] >= pLL_AF[0] + HingeGap
        # CheckUScase4 = perimeter[i,0]   >= pLL_AF[0] \
        #            and perimeter[i,0]   <  pLL_AF[0] + HingeGap \
        #            and perimeter[i-1,0] <  pLL_AF[0] + HingeGap
        CheckUScase5 = perimeter[i,0]   <  pLL_AF[0] \
                   and perimeter[i-1,0] >= pLL_AF[0] \
                   and perimeter[i-1,0] <  pLL_AF[0] + HingeGap
        CheckUScase6 = perimeter[i-1,0] <  pLL_AF[0]
        
        CheckLScase1 = perimeter[i,0]   <  pLL_AF[0]
        CheckLScase2 = perimeter[i,0]   >= pLL_AF[0] \
                   and cutMade == 0
        # CheckLScase3 = perimeter[i,0]   >= pLL_AF[0] \
        #            and cutMade == 1 \
        #            and perimeter[i,0]   <  pLL_AF[0] + cutWidth
        CheckLScase4 = perimeter[i,0]   >= pLL_AF[0] + cutWidth \
                   and cutMade == 1
        

        if CheckUS and CheckUScase1: 
            # US, i aft of ail hinge point, 1-1 mapping AF to AFwHinge
            j += 1
            perimeterWHinge[j,:] = perimeter[i,:]
            
        elif CheckUS and CheckUScase2:
            # US, i fore of both hinge points, i - 1 aft of both hinge points
            j += 1
            perimeterWHinge[j,0] = pLL_AF[0] + HingeGap
            perimeterWHinge[j,1] = np.interp(perimeterWHinge[j,0],AFGeom[:,0],AFGeom[:,1])
            jHingePoints[0] = j
            j += 1
            perimeterWHinge[j,0] = pLL_AF[0]
            perimeterWHinge[j,1] = np.interp(perimeterWHinge[j,0],AFGeom[:,0],AFGeom[:,1])
            jHingePoints[1] = j
            
        elif CheckUS and CheckUScase3:
            # US, i between hinge points, i-1 aft of ail hinge point 
            j += 1
            perimeterWHinge[j,0] = pLL_AF[0] + HingeGap
            perimeterWHinge[j,1] = np.interp(perimeterWHinge[j,0],AFGeom[:,0],AFGeom[:,1])
            jHingePoints[0] = j

        elif CheckUS and CheckUScase5:
            # US, i fore of wing hinge point, i-1 between hinge points
            j += 1
            perimeterWHinge[j,0] = pLL_AF[0]
            perimeterWHinge[j,1] = np.interp(perimeterWHinge[j,0],AFGeom[:,0],AFGeom[:,1])
            jHingePoints[1] = j
        elif CheckUS and CheckUScase6: 
            # US, i-1 fore of wing hinge point, 1-1 mapping AF to AFwHinge
            j += 1
            perimeterWHinge[j,:] = perimeter[i,:]
        elif CheckLS and CheckLScase1: 
            # LS, i fore of wing hinge point, 1-1 mapping AF to AFwHinge
            j += 1
            perimeterWHinge[j,:] = perimeter[i,:]
        elif CheckLS and CheckLScase2: 
            # LS, i first passed first hinge point, making hinge cut
            j += 1 
            perimeterWHinge[j,0] = pLL_AF[0]
            perimeterWHinge[j,1] = yLS_LL
            
            j += 1
            perimeterWHinge[j,0] = pLL_AF[0]
            perimeterWHinge[j,1] = yLS_LL + tauLL - (1-HingeOverlap)*extrWidth/chord
            jHingePoints[2] = j
            
            j += 1
            perimeterWHinge[j,0] = pLL_AF[0] + HingeGap
            yUSHP3 = np.interp(perimeterWHinge[j,0],AFGeom[:,0],AFGeom[:,1]) # yUS
            perimeterWHinge[j,1] = yUSHP3 - (1-HingeOverlap)*extrWidth/chord
            jHingePoints[3] = j
            
            j += 1
            perimeterWHinge[j,0] = pLL_AF[0] + cutWidth
            perimeterWHinge[j,1] = np.interp(perimeterWHinge[j,0],AFGeom[:,0],AFGeom[:,2]) # yLS
            
            i -= 1  # repeat ith point of perimeter
            cutMade = 1
            
        elif CheckLS and CheckLScase4: 
            # LS, i fore of wing hinge point, 1-1 mapping AF to AFwHinge
            j += 1
            perimeterWHinge[j,:] = perimeter[i,:]
        # else:
        #     print("Ail Point Error")
    
    perimeterWHinge = perimeterWHinge[0:j+1,:]
    nPerimWHinge = len(perimeterWHinge)
    
    # # Place LL at x-origin
    # perimeter[:,0] = perimeter[:,0] - pLL_AF[0] 
    # perimeter[:,1] = perimeter[:,1] - yLS_LL - tauLL*pLL_AF[1]
    
    # Scale and rotate to desired chord and angle of attack
    # perimeter = chord*perimeter
    # mRotAlpha = np.matrix([[np.cos(alpha),-np.sin(alpha)],
    #                     [np.sin(alpha), np.cos(alpha)]])
    # perimeter = perimeter@mRotAlpha
    
    # if iWing == 1: # Mirror about y-axis for Right Wing
    #     perimeter[:,0] = -perimeter[:,0]
    
    # # Translate to final position
    # perimeter[:,0] = perimeter[:,0] + pLL_bed[0]
    # perimeter[:,1] = perimeter[:,1] + pLL_bed[1]
    
    # Place LL at x-origin
    perimeterWHinge[:,0] = perimeterWHinge[:,0] - pLL_AF[0] 
    perimeterWHinge[:,1] = perimeterWHinge[:,1] - yLS_LL - tauLL*pLL_AF[1]
    
    perimeterWHinge = chord*perimeterWHinge
    mRotAlpha = np.matrix([[np.cos(alpha),-np.sin(alpha)],
                        [np.sin(alpha), np.cos(alpha)]])
    perimeterWHinge = perimeterWHinge@mRotAlpha
    
    if iWing == 1: # Mirror about y-axis for Right Wing
        perimeterWHinge[:,0] = -perimeterWHinge[:,0]
    
    # Translate to final position
    perimeterWHinge[:,0] = perimeterWHinge[:,0] + pLL_bed[0]
    perimeterWHinge[:,1] = perimeterWHinge[:,1] + pLL_bed[1]
    
    # if iWing == 1:
    #     xLE = float(max(perimeter[:,0]))
    # else:
    #     xLE = float(min(perimeter[:,0]))
    xMin = float(min(perimeterWHinge[:,0]))
    
    # Find US and LS Coords at Spar Locations -------------------------------
    AFPlaced = AFGeomAnalysis(perimeterWHinge,101)
    nSpars = len(xcStruct)
    xcStruct.sort() # Sort Structure LE -> TE
    xcStruct = np.flip(xcStruct) # Flip struture direction to TE -> LE
    if iWing == 1:
        xcStruct = -xcStruct + 1
        
    infill = np.zeros((nSpars,3))
    # for i in range(nSpars):
    #     infill[i][0] = chord*xcStruct[i] + xMin
    infill[:,0] = chord*xcStruct + xMin
    infill[:,1] = np.interp(infill[:,0],AFPlaced[:,0],AFPlaced[:,1]) - (1-infillOverlap)*extrWidth
    infill[:,2] = np.interp(infill[:,0],AFPlaced[:,0],AFPlaced[:,2]) + (1-infillOverlap)*extrWidth
    
    # breakpoint()
    
    # -------------------- Perimeter GCode Composition ---------------------- #
    # Initialize perimeter gcode
    # breakpoint()           
    perimeterHingeStart = np.concatenate((perimeterWHinge[int(jHingePoints[3]):nPerimWHinge,:], perimeterWHinge[0:int(jHingePoints[3]),:]), axis=0)
    perimeterR = np.flipud(perimeterHingeStart) # reverse perimeter to out LS then US
    # perimeterR = perimeterHingeStart # Uncomment for Flap then Main element
    nPerimHS = len(perimeterR)
    perimeter_command = "\n\nM117 L%.0f/%.0f;Z%.1f"% (iLyrs,nLyrs,Z)
    # breakpoint()
    # Move to Perimeter start
    cmndG = "\nG0"
    cmndF = " F%.2f" % (60*speedTravel)
    cmndX = " X%.4f" % (perimeterR[0,0]+iWing*flangeLength*extrWidth)
    cmndY = " Y%.4f" % (perimeterR[0,1])
    cmndZ = " Z%.4f" % (Z+Zhop)
    perimeter_command = perimeter_command + cmndG
    perimeter_command = perimeter_command + cmndF
    perimeter_command = perimeter_command + cmndX
    perimeter_command = perimeter_command + cmndY
    perimeter_command = perimeter_command + cmndZ
    perimeter_command = perimeter_command + "  ;Move to Perimeter Start"
    
    # Set Z Height
    cmndG = "\nG0"
    cmndF = " F%.2f" % (60*speedTravel)
    cmndX = " X%.4f" % (perimeterR[0,0]+iWing*flangeLength*extrWidth)
    cmndY = " Y%.4f" % (perimeterR[0,1])
    cmndZ = " Z%.4f" % (Z)
    perimeter_command = perimeter_command + cmndG
    perimeter_command = perimeter_command + cmndF
    perimeter_command = perimeter_command + cmndX
    perimeter_command = perimeter_command + cmndY
    perimeter_command = perimeter_command + cmndZ
    perimeter_command = perimeter_command + "  ;Set Z Height"
    
    # Prime Nozzle
    if iLyrs == 0:
        Eprime = extrWidth*layer0Height
        dE = FrstLyrMult*FlowMult*Eprime*flangeLength/FilaDiameter**2
    else:
        Eprime = extrWidth*layerHeight
        dE = FlowMult*Eprime*flangeLength/FilaDiameter**2
    E = E + Ehop + dE
    cmndG = "\nG1"
    cmndF = " F%.2f" % (60*speedExtrusion)
    cmndE = " E%.2f" % (E)
    cmndX = " X%.4f" % (perimeterR[0,0])
    cmndY = " Y%.4f" % (perimeterR[0,1])
    cmndZ = " Z%.4f" % (Z)
    perimeter_command = perimeter_command + cmndG
    perimeter_command = perimeter_command + cmndF
    perimeter_command = perimeter_command + cmndE
    perimeter_command = perimeter_command + cmndX
    perimeter_command = perimeter_command + cmndY
    perimeter_command = perimeter_command + cmndZ
    perimeter_command = perimeter_command + "  ;Prime Nozzle"
    
    for i1 in range(1,nPerimHS):
        # atHinge = i1 == jHingePoints[0] + nPerimHS - jHingePoints[3]
        atHinge = i1 == - jHingePoints[0] + jHingePoints[3] - 1
        if atHinge:
            cmndG = "\nG0"
        else:
            cmndG = "\nG1"
        
        nzlXo = perimeterR[i1-1,0]
        nzlYo = perimeterR[i1-1,1]
        nzlX  = perimeterR[i1  ,0]
        nzlY  = perimeterR[i1  ,1]
        nzlTrvl = np.sqrt((nzlX - nzlXo)**2 + (nzlY - nzlYo)**2)
        cmndX = " X%.4f" % (nzlX)
        cmndY = " Y%.4f" % (nzlY)
        
        if atHinge == 0:
            if iLyrs == 0:
                Eprime = extrWidth*layer0Height
                dE = FrstLyrMult*FlowMult*Eprime*nzlTrvl/FilaDiameter**2
            else:
                Eprime = extrWidth*layerHeight
                dE = FlowMult*Eprime*nzlTrvl/FilaDiameter**2
                
            E = E + dE
            cmndE = " E%.4f" % (E)
        
        perimeter_command = perimeter_command + cmndG
        perimeter_command = perimeter_command + cmndF
        if atHinge == 0:
            perimeter_command = perimeter_command + cmndE
        perimeter_command = perimeter_command + cmndX
        perimeter_command = perimeter_command + cmndY
        if atHinge:
            perimeter_command = perimeter_command + "  ;Hinge"
        else:
            perimeter_command = perimeter_command + "  ;Perimeter"
    
    # cmndG = "\nG1"
    
    # nzlXo = perimeterR[-1,0]
    # nzlYo = perimeterR[-1,1]
    # nzlX  = perimeterR[ 0,0]
    # nzlY  = perimeterR[ 0,1]
    
    # nzlTrvl = np.sqrt((nzlX - nzlXo)**2 + (nzlY - nzlYo)**2)
    # if iLyrs == 0:
    #     Eprime = extrWidth*layer0Height
    #     dE = FrstLyrMult*FlowMult*Eprime*nzlTrvl/FilaDiameter**2
    # else:
    #     Eprime = extrWidth*layerHeight
    #     dE = FlowMult*Eprime*nzlTrvl/FilaDiameter**2
        
    # E = E + dE
    # cmndE = " E%.4f" % (E)
    
    # cmndX = " X%.4f" % (perimeterR[0,0])
    # cmndY = " Y%.4f" % (perimeterR[0,1])
    # perimeter_command = perimeter_command + cmndG
    # perimeter_command = perimeter_command + cmndF
    # perimeter_command = perimeter_command + cmndE
    # perimeter_command = perimeter_command + cmndX
    # perimeter_command = perimeter_command + cmndY
    # perimeter_command = perimeter_command + cmndZ
    # perimeter_command = perimeter_command + "  ;Perimeter End"
            
    fo.write(perimeter_command)
    
    # ---------------- Internal Structure GCode Composition ----------------- #
    sparDirection = 0 # 0 = always top to bottom, 1 = alternate directions
    
    infill_command = "\n"
    
    # Lift Nozzle with retraction
    E = E - Ehop
    cmndG = "\nG1"
    cmndF = " F%.2f" % (60*speedExtrusion)
    cmndE = " E%.2f" % (E)
    cmndZ = " Z%.4f" % (Z+Zhop)
    infill_command = infill_command + cmndG
    infill_command = infill_command + cmndF
    infill_command = infill_command + cmndE
    infill_command = infill_command + cmndX
    infill_command = infill_command + cmndY
    infill_command = infill_command + cmndZ
    infill_command = infill_command + "  ;Retraction/Zhop"
    
    for i in range(nSpars):
        # Travel to Spar Location
        cmndG = "\nG0"
        cmndF = " F%.2f" % (60*speedTravel)
        cmndX = " X%.4f" % (infill[i][0] - iWing*flangeLength*extrWidth)
        if sparDirection: # alternate directions
            if np.mod(i,2):
                cmndY = " Y%.4f" % (infill[i][1])
            else:
                cmndY = " Y%.4f" % (infill[i][2])
        else: # always upper to lower
            cmndY = " Y%.4f" % (infill[i][1])
        infill_command = infill_command + cmndG
        infill_command = infill_command + cmndF
        infill_command = infill_command + cmndX
        infill_command = infill_command + cmndY
        infill_command = infill_command + cmndZ
        infill_command = infill_command + "  ;Travel to Spar"
        
        # Set Z Height
        cmndZ = " Z%.4f" % (Z)
        infill_command = infill_command + cmndG
        infill_command = infill_command + cmndF
        infill_command = infill_command + cmndX
        infill_command = infill_command + cmndY
        infill_command = infill_command + cmndZ
        infill_command = infill_command + "  ;Set Z Height"
        
        # Prime Nozzle
        if iLyrs == 0:
            Eprime = extrWidth*layer0Height
            dE = FrstLyrMult*FlowMult*Eprime*flangeLength/FilaDiameter**2
        else:
            Eprime = extrWidth*layerHeight
            dE = FlowMult*Eprime*flangeLength/FilaDiameter**2
        E = E + Ehop + dE
        cmndG = "\nG1"
        cmndF = " F%.2f" % (60*speedExtrusion)
        cmndE = " E%.2f" % (E)
        cmndX = " X%.4f" % (infill[i][0])
        infill_command = infill_command + cmndG
        infill_command = infill_command + cmndF
        infill_command = infill_command + cmndE
        infill_command = infill_command + cmndX
        infill_command = infill_command + cmndY
        infill_command = infill_command + cmndZ
        infill_command = infill_command + "  ;Prime Nozzle"
        
        # Spar Command
        nzlTrvl = np.sqrt((infill[i][1] - infill[i][2])**2)
        if iLyrs == 0:
            Eprime = extrWidth*layer0Height
            dE = FrstLyrMult*FlowMult*Eprime*nzlTrvl/FilaDiameter**2
        else:
            Eprime = extrWidth*layerHeight
            dE = FlowMult*Eprime*nzlTrvl/FilaDiameter**2
        E = E + dE
        cmndF = " F%.2f" % (60*speedInfExtrusion)
        cmndE = " E%.2f" % (E)
        cmndX = " X%.4f" % (infill[i][0])
        if sparDirection: # alternate directions
            if np.mod(i,2):
                cmndY = " Y%.4f" % (infill[i][2])
            else:
                cmndY = " Y%.4f" % (infill[i][1])
        else: # always upper to lower
            cmndY = " Y%.4f" % (infill[i][2])
        infill_command = infill_command + cmndG
        infill_command = infill_command + cmndF
        infill_command = infill_command + cmndE
        infill_command = infill_command + cmndX
        infill_command = infill_command + cmndY
        infill_command = infill_command + cmndZ
        infill_command = infill_command + "  ;Spar"
        
        # Lift Nozzle with retraction
        E = E - Ehop
        cmndG = "\nG1"
        cmndE = " E%.2f" % (E)
        cmndZ = " Z%.4f" % (Z+Zhop)
        infill_command = infill_command + cmndG
        infill_command = infill_command + cmndF
        infill_command = infill_command + cmndE
        infill_command = infill_command + cmndX
        infill_command = infill_command + cmndY
        infill_command = infill_command + cmndZ
        infill_command = infill_command + "  ;Retraction/Zhop"
    fo.write(infill_command)
    return E

def RibWAil(fo,iLyrs,nLyrs,AFcoords,chord,alpha,pLL_AF,pLL_bed,extrWidth,speedExtrusion,speedInfExtrusion,speedTravel,layer0Height,layerHeight,E,Z,Ehop,EhopP,Zhop,FlowMult,FrstLyrMult,SldFillMult,FilaDiameter,infillOverlap,iWing,HingeGap,HingeOverlap):
    # breakpoint()
    import numpy as np
    from CodeFiles.Aerodynamics.Airfoils.operations import AFGeomAnalysis
    chord = float(chord) # ensuring chord is acutally a float
    nperim = len(AFcoords)  # Define perimeter length
    
    # Define initial perimeter ------------------
    perimeter = np.array(AFcoords)
    
    # TE thickness
    
    # rediscretize AF from Selig to [x,yUS,yLS,t,h]
    AFGeom = AFGeomAnalysis(perimeter,201)
    tauLL  = np.interp(pLL_AF[0],AFGeom[:,0],AFGeom[:,3]) # tau @LL
    yLS_LL = np.interp(pLL_AF[0],AFGeom[:,0],AFGeom[:,2]) # yLS @LL
    
    # generate airfoil with hinge
    perimeterWHinge = np.zeros((2*nperim,2))
    jHingePoints = np.zeros(4)
    HingeGap = HingeGap*extrWidth/chord
    cutWidth   = HingeGap + tauLL
    cutMade = 0
    CheckUS = 1
    CheckLS = 0
    i = -1 # initialize counter for base perimeter
    j = -1 # initialize counter for perimeter with hinge
    while i < nperim - 1:
        i += 1
        if i >2:
            CheckUS = perimeter[i,0] <  perimeter[i-1,0]
            CheckLS = perimeter[i,0] >= perimeter[i-1,0]
        
        CheckUScase1 = perimeter[i,0]   >= pLL_AF[0] + HingeGap
        CheckUScase2 = perimeter[i,0]   <  pLL_AF[0] \
                   and perimeter[i-1,0] >= pLL_AF[0] + HingeGap
        CheckUScase3 = perimeter[i,0]   >= pLL_AF[0] \
                   and perimeter[i,0]   <  pLL_AF[0] + HingeGap \
                   and perimeter[i-1,0] >= pLL_AF[0] + HingeGap
        CheckUScase3 = perimeter[i,0]   >= pLL_AF[0] \
                   and perimeter[i,0]   <  pLL_AF[0] + HingeGap \
                   and perimeter[i-1,0] >= pLL_AF[0] + HingeGap
        # CheckUScase4 = perimeter[i,0]   >= pLL_AF[0] \
        #            and perimeter[i,0]   <  pLL_AF[0] + HingeGap \
        #            and perimeter[i-1,0] <  pLL_AF[0] + HingeGap
        CheckUScase5 = perimeter[i,0]   <  pLL_AF[0] \
                   and perimeter[i-1,0] >= pLL_AF[0] \
                   and perimeter[i-1,0] <  pLL_AF[0] + HingeGap
        CheckUScase6 = perimeter[i-1,0] <  pLL_AF[0]
        
        CheckLScase1 = perimeter[i,0]   <  pLL_AF[0]
        CheckLScase2 = perimeter[i,0]   >= pLL_AF[0] \
                   and cutMade == 0
        # CheckLScase3 = perimeter[i,0]   >= pLL_AF[0] \
        #            and cutMade == 1 \
        #            and perimeter[i,0]   <  pLL_AF[0] + cutWidth
        CheckLScase4 = perimeter[i,0]   >= pLL_AF[0] + cutWidth \
                   and cutMade == 1
        

        if CheckUS and CheckUScase1: 
            # US, i aft of ail hinge point, 1-1 mapping AF to AFwHinge
            j += 1
            perimeterWHinge[j,:] = perimeter[i,:]
            
        elif CheckUS and CheckUScase2:
            # US, i fore of both hinge points, i - 1 aft of both hinge points
            j += 1
            perimeterWHinge[j,0] = pLL_AF[0] + HingeGap
            perimeterWHinge[j,1] = np.interp(perimeterWHinge[j,0],AFGeom[:,0],AFGeom[:,1])
            jHingePoints[0] = j
            j += 1
            perimeterWHinge[j,0] = pLL_AF[0]
            perimeterWHinge[j,1] = np.interp(perimeterWHinge[j,0],AFGeom[:,0],AFGeom[:,1])
            jHingePoints[1] = j
            
        elif CheckUS and CheckUScase3:
            # US, i between hinge points, i-1 aft of ail hinge point 
            j += 1
            perimeterWHinge[j,0] = pLL_AF[0] + HingeGap
            perimeterWHinge[j,1] = np.interp(perimeterWHinge[j,0],AFGeom[:,0],AFGeom[:,1])
            jHingePoints[0] = j

        elif CheckUS and CheckUScase5:
            # US, i fore of wing hinge point, i-1 between hinge points
            j += 1
            perimeterWHinge[j,0] = pLL_AF[0]
            perimeterWHinge[j,1] = np.interp(perimeterWHinge[j,0],AFGeom[:,0],AFGeom[:,1])
            jHingePoints[1] = j
        elif CheckUS and CheckUScase6: 
            # US, i-1 fore of wing hinge point, 1-1 mapping AF to AFwHinge
            j += 1
            perimeterWHinge[j,:] = perimeter[i,:]
        elif CheckLS and CheckLScase1: 
            # LS, i fore of wing hinge point, 1-1 mapping AF to AFwHinge
            j += 1
            perimeterWHinge[j,:] = perimeter[i,:]
        elif CheckLS and CheckLScase2: 
            # LS, i first passed first hinge point, making hinge cut
            j += 1 
            perimeterWHinge[j,0] = pLL_AF[0]
            perimeterWHinge[j,1] = yLS_LL
            
            j += 1
            perimeterWHinge[j,0] = pLL_AF[0]
            perimeterWHinge[j,1] = yLS_LL + tauLL - HingeOverlap*extrWidth/chord
            jHingePoints[2] = j
            
            j += 1
            perimeterWHinge[j,0] = pLL_AF[0] + HingeGap
            yUSHP3 = np.interp(perimeterWHinge[j,0],AFGeom[:,0],AFGeom[:,1]) # yUS
            perimeterWHinge[j,1] = yUSHP3 - HingeOverlap*extrWidth/chord
            jHingePoints[3] = j
            
            j += 1
            perimeterWHinge[j,0] = pLL_AF[0] + cutWidth
            perimeterWHinge[j,1] = np.interp(perimeterWHinge[j,0],AFGeom[:,0],AFGeom[:,2]) # yLS
            
            i -= 1  # repeat ith point of perimeter
            cutMade = 1
            
        elif CheckLS and CheckLScase4: 
            # LS, i fore of wing hinge point, 1-1 mapping AF to AFwHinge
            j += 1
            perimeterWHinge[j,:] = perimeter[i,:]
        # else:
        #     print("Ail Point Error")
    
    perimeterWHinge = perimeterWHinge[0:j+1,:]
    nPerimWHinge = len(perimeterWHinge)
    
    # # Place LL at x-origin
    # perimeter[:,0] = perimeter[:,0] - pLL_AF[0] 
    # perimeter[:,1] = perimeter[:,1] - yLS_LL - tauLL*pLL_AF[1]
    
    # Scale and rotate to desired chord and angle of attack
    # perimeter = chord*perimeter
    # mRotAlpha = np.matrix([[np.cos(alpha),-np.sin(alpha)],
    #                     [np.sin(alpha), np.cos(alpha)]])
    # perimeter = perimeter@mRotAlpha
    
    # if iWing == 1: # Mirror about y-axis for Right Wing
    #     perimeter[:,0] = -perimeter[:,0]
    
    # # Translate to final position
    # perimeter[:,0] = perimeter[:,0] + pLL_bed[0]
    # perimeter[:,1] = perimeter[:,1] + pLL_bed[1]
    
    # Place LL at x-origin
    perimeterWHinge[:,0] = perimeterWHinge[:,0] - pLL_AF[0] 
    perimeterWHinge[:,1] = perimeterWHinge[:,1] - yLS_LL - tauLL*pLL_AF[1]
    
    perimeterWHinge = chord*perimeterWHinge
    mRotAlpha = np.matrix([[np.cos(alpha),-np.sin(alpha)],
                        [np.sin(alpha), np.cos(alpha)]])
    perimeterWHinge = perimeterWHinge@mRotAlpha
    
    if iWing == 1: # Mirror about y-axis for Right Wing
        perimeterWHinge[:,0] = -perimeterWHinge[:,0]
    
    # Translate to final position
    perimeterWHinge[:,0] = perimeterWHinge[:,0] + pLL_bed[0]
    perimeterWHinge[:,1] = perimeterWHinge[:,1] + pLL_bed[1]
    
    # if iWing == 1:
    #     xLE = float(max(perimeter[:,0]))
    # else:
    #     xLE = float(min(perimeter[:,0]))
    xMin = float(min(perimeterWHinge[:,0]))
    xMax = float(max(perimeterWHinge[:,0]))
    
    # # Find US and LS Coords at Spar Locations -------------------------------
    # AFPlaced = AFGeomAnalysis(perimeterWHinge,101)
    # nSpars = len(xcStruct)
    # xcStruct.sort() # Sort Structure LE -> TE
    # xcStruct = np.flip(xcStruct) # Flip struture direction to TE -> LE
    # if iWing == 1:
    #     xcStruct = -xcStruct + 1
        
    # infill = np.zeros((nSpars,3))
    # # for i in range(nSpars):
    # #     infill[i][0] = chord*xcStruct[i] + xMin
    # infill[:,0] = chord*xcStruct + xMin
    # infill[:,1] = np.interp(infill[:,0],AFPlaced[:,0],AFPlaced[:,1]) - infillOverlap*extrWidth
    # infill[:,2] = np.interp(infill[:,0],AFPlaced[:,0],AFPlaced[:,2]) + infillOverlap*extrWidth
    
    # # breakpoint()
    
    # Solid Fill Calculations -----------------------------------------------
    AFPlaced = AFGeomAnalysis(perimeterWHinge,101)
    solidFillLEBuffer = 5
    solidFillTEBuffer = 30
    if iWing == 1:
        xLE = xMax
        xTE = xMin
    else:
        xLE = xMin
        xTE = xMax
        
    infill = np.zeros((1,3))
    infillStart = xTE + iWing*extrWidth*solidFillTEBuffer
    infill[0,0] = infillStart
    infill[0,1] = np.interp(infill[:,0],AFPlaced[:,0],AFPlaced[:,1]) - infillOverlap*extrWidth
    infill[0,2] = np.interp(infill[:,0],AFPlaced[:,0],AFPlaced[:,2]) + infillOverlap*extrWidth
    i = 0
    while 1:
        i += 1
        infill = np.append(infill,[[0,0,0]],axis = 0)
        infill[i,0] = infill[i-1,0] + iWing*extrWidth
        infill[i,1] = np.interp(infill[i,0],AFPlaced[:,0],AFPlaced[:,1]) - infillOverlap*extrWidth
        infill[i,2] = np.interp(infill[i,0],AFPlaced[:,0],AFPlaced[:,2]) + infillOverlap*extrWidth
        if iWing*(xLE - iWing*solidFillLEBuffer*extrWidth - infill[i,0]) < 0:
            break
            
    nFill = len(infill[:,0])
    
    # -------------------- Perimeter GCode Composition ---------------------- #
    # Initialize perimeter gcode
    # breakpoint()           
    perimeterHingeStart = np.concatenate((perimeterWHinge[int(jHingePoints[3]):nPerimWHinge,:], perimeterWHinge[0:int(jHingePoints[3]),:]), axis=0)
    perimeterR = np.flipud(perimeterHingeStart) # reverse perimeter to out LS then US
    nPerimHS = len(perimeterR)
    perimeter_command = "\n\nM117 L%.0f/%.0f;Z%.1f"% (iLyrs,nLyrs,Z)
    # breakpoint()
    # Move to Perimeter start
    cmndG = "\nG0"
    cmndF = " F%.2f" % (60*speedTravel)
    cmndX = " X%.4f" % (perimeterR[0,0])
    cmndY = " Y%.4f" % (perimeterR[0,1])
    cmndZ = " Z%.4f" % (Z+Zhop)
    perimeter_command = perimeter_command + cmndG
    perimeter_command = perimeter_command + cmndF
    perimeter_command = perimeter_command + cmndX
    perimeter_command = perimeter_command + cmndY
    perimeter_command = perimeter_command + cmndZ
    perimeter_command = perimeter_command + "  ;Move to Perimeter Start"
    
    # Set Z Height
    cmndG = "\nG0"
    cmndF = " F%.2f" % (60*speedTravel)
    cmndX = " X%.4f" % (perimeterR[0,0])
    cmndY = " Y%.4f" % (perimeterR[0,1])
    cmndZ = " Z%.4f" % (Z)
    perimeter_command = perimeter_command + cmndG
    perimeter_command = perimeter_command + cmndF
    perimeter_command = perimeter_command + cmndX
    perimeter_command = perimeter_command + cmndY
    perimeter_command = perimeter_command + cmndZ
    perimeter_command = perimeter_command + "  ;Set Z Height"
    
    # Prime Nozzle
    E = E + Ehop + EhopP
    cmndG = "\nG1"
    cmndF = " F%.2f" % (60*speedExtrusion)
    cmndE = " E%.2f" % (E)
    cmndZ = " Z%.4f" % (Z)
    perimeter_command = perimeter_command + cmndG
    perimeter_command = perimeter_command + cmndF
    perimeter_command = perimeter_command + cmndE
    perimeter_command = perimeter_command + cmndX
    perimeter_command = perimeter_command + cmndY
    perimeter_command = perimeter_command + cmndZ
    perimeter_command = perimeter_command + "  ;Prime Nozzle"
    
    for i1 in range(1,nPerimHS):
        # atHinge = i1 == jHingePoints[0] + nPerimHS - jHingePoints[3]
        atHinge = i1 == - jHingePoints[0] + jHingePoints[3] - 1
        if atHinge:
            cmndG = "\nG0"
        else:
            cmndG = "\nG1"
        
        nzlXo = perimeterR[i1-1,0]
        nzlYo = perimeterR[i1-1,1]
        nzlX  = perimeterR[i1  ,0]
        nzlY  = perimeterR[i1  ,1]
        nzlTrvl = np.sqrt((nzlX - nzlXo)**2 + (nzlY - nzlYo)**2)
        cmndX = " X%.4f" % (nzlX)
        cmndY = " Y%.4f" % (nzlY)
        
        if atHinge == 0:
            if iLyrs == 0:
                Eprime = extrWidth*layer0Height
                dE = FrstLyrMult*FlowMult*Eprime*nzlTrvl/FilaDiameter**2
            else:
                Eprime = extrWidth*layerHeight
                dE = FlowMult*Eprime*nzlTrvl/FilaDiameter**2
                
            E = E + dE
            cmndE = " E%.4f" % (E)
        
        perimeter_command = perimeter_command + cmndG
        perimeter_command = perimeter_command + cmndF
        if atHinge == 0:
            perimeter_command = perimeter_command + cmndE
        perimeter_command = perimeter_command + cmndX
        perimeter_command = perimeter_command + cmndY
        if atHinge:
            perimeter_command = perimeter_command + "  ;Hinge"
        else:
            perimeter_command = perimeter_command + "  ;Perimeter"
    
            
    fo.write(perimeter_command)
    
    # ---------------- Rib (Solid Fill) GCode Composition ----------------- #
    infill_command = "\n"
    
    # # Lift Nozzle with retraction
    # E = E - Ehop
    # cmndG = "\nG1"
    # cmndF = " F%.2f" % (60*speedExtrusion)
    # cmndE = " E%.2f" % (E)
    # cmndZ = " Z%.4f" % (Z+Zhop)
    # infill_command = infill_command + cmndG
    # infill_command = infill_command + cmndF
    # infill_command = infill_command + cmndE
    # infill_command = infill_command + cmndX
    # infill_command = infill_command + cmndY
    # infill_command = infill_command + cmndZ
    # infill_command = infill_command + "  ;Retraction/Zhop"
    # breakpoint()
    for i in range(nFill):
        # Travel to Spar Location
        cmndG = "\nG0"
        cmndF = " F%.2f" % (60*speedTravel)
        cmndX = " X%.4f" % (infill[i][0])
        if np.mod(i,2):
            cmndY = " Y%.4f" % (infill[i][1])
        else:
            cmndY = " Y%.4f" % (infill[i][2])
        infill_command = infill_command + cmndG
        infill_command = infill_command + cmndF
        infill_command = infill_command + cmndX
        infill_command = infill_command + cmndY
        infill_command = infill_command + cmndZ
        infill_command = infill_command + "  ;Travel to Fill"
        
        # # Set Z Height
        # cmndZ = " Z%.4f" % (Z)
        # infill_command = infill_command + cmndG
        # infill_command = infill_command + cmndF
        # infill_command = infill_command + cmndX
        # infill_command = infill_command + cmndY
        # infill_command = infill_command + cmndZ
        # infill_command = infill_command + "  ;Set Z Height"
        
        # # Prime Nozzle
        # E = E + Ehop
        # cmndG = "\nG1"
        # cmndF = " F%.2f" % (60*speedExtrusion)
        # cmndE = " E%.2f" % (E)
        # infill_command = infill_command + cmndG
        # infill_command = infill_command + cmndF
        # infill_command = infill_command + cmndE
        # infill_command = infill_command + cmndX
        # infill_command = infill_command + cmndY
        # infill_command = infill_command + cmndZ
        # infill_command = infill_command + "  ;Prime Nozzle"
        
        # Spar Command
        nzlTrvl = np.sqrt((infill[i][1] - infill[i][2])**2)
        Eprime = extrWidth*layer0Height
        dE = FrstLyrMult*FlowMult*Eprime*nzlTrvl/FilaDiameter**2
        E = E + dE
        cmndF = " F%.2f" % (60*speedInfExtrusion)
        cmndE = " E%.2f" % (E)
        cmndX = " X%.4f" % (infill[i][0])
        if np.mod(i,2):
            cmndY = " Y%.4f" % (infill[i][2])
        else:
            cmndY = " Y%.4f" % (infill[i][1])
        infill_command = infill_command + cmndG
        infill_command = infill_command + cmndF
        infill_command = infill_command + cmndE
        infill_command = infill_command + cmndX
        infill_command = infill_command + cmndY
        infill_command = infill_command + cmndZ
        infill_command = infill_command + "  ;Rib"
        
    # Lift Nozzle with retraction
    E = E - Ehop
    cmndG = "\nG1"
    cmndE = " E%.2f" % (E)
    cmndZ = " Z%.4f" % (Z+Zhop)
    infill_command = infill_command + cmndG
    infill_command = infill_command + cmndF
    infill_command = infill_command + cmndE
    infill_command = infill_command + cmndX
    infill_command = infill_command + cmndY
    infill_command = infill_command + cmndZ
    infill_command = infill_command + "  ;Retraction/Zhop"
    fo.write(infill_command)
    return E
