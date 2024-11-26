# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 13:33:55 2020

@author: justi
"""

# ------------------------- Brim Generation -------------------------- #

nperim = len(airfoil[:,1])
perimeter0 = np.array(np.zeros((nperim,2)))
perimeter0[:,0] = np.array(c_airfoil*airfoil[:,0] + origin[0])
perimeter0[:,1] = np.array(c_airfoil*airfoil[:,1] + origin[1])
 
Skirts = np.array(np.zeros((nperim,2,skirtnum)))
for i2 in list(range(skirtnum)):
    for i in list(range(nperim)): 
        if i2 == 0:
            if i == 0:
                P0 = np.array([[perimeter0[nperim-1,0]],[perimeter0[nperim-1,1]]])
                P1 = np.array([[perimeter0[     0 ,0]],[perimeter0[     0 ,1]]])
                P2 = np.array([[perimeter0[     1 ,0]],[perimeter0[     1 ,1]]])
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

        Phi = 0*np.pi/180
        MrotPhi = np.matrix([[np.cos(Phi),-np.sin(Phi)],
                           [np.sin(Phi), np.cos(Phi)]])

        P0 = np.matmul(MrotPhi,P0)
        P1 = np.matmul(MrotPhi,P1)
        P2 = np.matmul(MrotPhi,P2)

        Mrot90 = np.matrix([[np.cos(-np.pi/2),-np.sin(-np.pi/2)],
                            [np.sin(-np.pi/2), np.cos(-np.pi/2)]])
        t0 = P1 - P0
        t0 = t0/np.linalg.norm(t0)
        t2 = P2 - P1
        t2 = t2/np.linalg.norm(t2)

        n0 = np.matmul(Mrot90,t0)
        n2 = np.matmul(Mrot90,t2)

        n1 = (n0 + n2)
        n1 = n1/np.linalg.norm(n1)

        theta0 = np.arctan2(t0[1,0],t0[0,0])
        MrotTheta0 = np.matrix([[np.cos(-theta0),-np.sin(-theta0)],
                                [np.sin(-theta0), np.cos(-theta0)]])
        dt = np.matmul(MrotTheta0,t2)

        theta = np.arctan2(dt[1,0],dt[0,0])
        if i2 == 0:
            eps_mag = (0.5*extrWidth + skirtOffset)/(np.cos(theta/2))
        else:
            eps_mag =     (extrWidth)/(np.cos(theta/2))
        eps = eps_mag*n1
        P1offset = P1 + eps

        Skirts[i,0,i2] = P1offset[0,0]
        Skirts[i,1,i2] = P1offset[1,0]

# ------------------------ Brim GCode Composition -------------------- #
brim_command = "\n\nM117Skirt\n"
cmndG = "G0"
cmndZ = " Z%.4f" % (Z)
cmndF = " F%.2f" % (60*speedExtrusion)
brim_command = brim_command + cmndG
brim_command = brim_command + cmndZ
brim_command = brim_command + cmndF
brim_command = brim_command + "  ;Set Z"

for i2 in list(range(skirtnum)):
    brim_command = brim_command + "\n"
    cmndG = "G0"
    cmndX = " X%.4f" % (Skirts[0,0,skirtnum - 1 - i2])
    cmndY = " Y%.4f" % (Skirts[0,1,skirtnum - 1 - i2])
    brim_command = brim_command + cmndG
    brim_command = brim_command + cmndF
    brim_command = brim_command + cmndX
    brim_command = brim_command + cmndY
    brim_command = brim_command + "  ;Brim Loop Start"


    for i1 in list(range(1,nperim)):
        brim_command = brim_command + "\n"
        cmndG = "G1"

        nzlXo = Skirts[i1-1,0,skirtnum - 1 - i2]
        nzlYo = Skirts[i1-1,1,skirtnum - 1 - i2]
        nzlX = Skirts[i1,0,skirtnum - 1 - i2]
        nzlY = Skirts[i1,1,skirtnum - 1 - i2]
        nzlTrvl = np.sqrt((nzlX - nzlXo)**2 + (nzlY - nzlYo)**2)
        cmndX = " X%.4f" % (nzlX)
        cmndY = " Y%.4f" % (nzlY)

        if i2 != 0:
            Eprime = extrWidth*layer0Height
        else:
            Eprime = extrWidth*layerHeight
        dE = SldFillMult*FlowMult*Eprime*nzlTrvl/FilaDiameter**2
        E = E + dE
        cmndE = " E%.4f" % (E)

        brim_command = brim_command + cmndG
        brim_command = brim_command + cmndF
        brim_command = brim_command + cmndE
        brim_command = brim_command + cmndX
        brim_command = brim_command + cmndY
        brim_command = brim_command + "  ;Brim"

    brim_command = brim_command + "\n"
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
    brim_command = brim_command + cmndG
    brim_command = brim_command + cmndF
    brim_command = brim_command + cmndE
    brim_command = brim_command + cmndX
    brim_command = brim_command + cmndY
    brim_command = brim_command + "  ;Brim"
        
# ------------------------ Perimeters ------------------------- #
hull = np.array(np.zeros((len(airfoil[:,0]),2)))
perimeters = np.array(np.zeros((len(airfoil[:,0]),2,n_lyrs)))
perimeters_command = [None]*n_lyrs

"""
the i2 loop iterates through every layer, offsets the toolpath so
that the part is circumscribed by the airfoil, and composes the 
perimeter gcode"""
for i2 in list(range(n_lyrs)):
    
    hull[:,0] = c_airfoil*airfoil[:,0] + origin[0]
    hull[:,1] = c_airfoil*airfoil[:,1] + origin[1]
    nhull = len(hull[:,0])
    
    
    for i in list(range(nhull)):    
        if i == 0:
            P0 = np.array([[hull[nhull-1,0]],[hull[nhull-1,1]]])
            P1 = np.array([[hull[     0 ,0]],[hull[     0 ,1]]])
            P2 = np.array([[hull[     1 ,0]],[hull[     1 ,1]]])
        elif i == nhull - 1:
            P0 = np.array([[hull[i-1,0]],[hull[i-1,1]]])
            P1 = np.array([[hull[ i ,0]],[hull[ i ,1]]])
            P2 = np.array([[hull[ 0 ,0]],[hull[ 0 ,1]]])
        else:
            P0 = np.array([[hull[i-1,0]],[hull[i-1,1]]])
            P1 = np.array([[hull[ i ,0]],[hull[ i ,1]]])
            P2 = np.array([[hull[i+1,0]],[hull[i+1,1]]])

        Phi = 0*np.pi/180
        MrotPhi = np.matrix([[np.cos(Phi),-np.sin(Phi)],
                           [np.sin(Phi), np.cos(Phi)]])

        P0 = np.matmul(MrotPhi,P0)
        P1 = np.matmul(MrotPhi,P1)
        P2 = np.matmul(MrotPhi,P2)

        Mrot90 = np.matrix([[np.cos(-np.pi/2),-np.sin(-np.pi/2)],
                           [np.sin(-np.pi/2), np.cos(-np.pi/2)]])
        t0 = P1 - P0
        t0 = t0/np.linalg.norm(t0)
        t2 = P2 - P1
        t2 = t2/np.linalg.norm(t2)

        n0 = np.matmul(Mrot90,t0)
        n2 = np.matmul(Mrot90,t2)

        n1 = (n0 + n2)
        n1 = n1/np.linalg.norm(n1)

        theta0 = np.arctan2(t0[1,0],t0[0,0])
        MrotTheta0 = np.matrix([[np.cos(-theta0),-np.sin(-theta0)],
                                [np.sin(-theta0), np.cos(-theta0)]])
        dt = np.matmul(MrotTheta0,t2)

        theta = np.arctan2(dt[1,0],dt[0,0])
        eps_mag = -(0.5*extrWidth)/(np.cos(theta/2))
        eps = eps_mag*n1
        P1offset = P1 + eps

        perimeters[i,0,i2] = P1offset[0,0]
        perimeters[i,1,i2] = P1offset[1,0]
    
    if (i2 != 0):
        Z = Z + layerHeight
        
    perimeters_command[i2] = "\n\nM117 L%.0f/%.0f,Z%.1f"% (i2,n_lyrs,Z)
    
    cmndG = "\nG0"
    cmndX = " X%.4f" % (perimeters[0,0,i2])
    cmndY = " Y%.4f" % (perimeters[0,1,i2])
    perimeters_command[i2] = perimeters_command[i2] + cmndG
    perimeters_command[i2] = perimeters_command[i2] + cmndF
    perimeters_command[i2] = perimeters_command[i2] + cmndX
    perimeters_command[i2] = perimeters_command[i2] + cmndY
    perimeters_command[i2] = perimeters_command[i2] + "  ;Perimeter Start"
    
    cmndG = "\nG0"
    cmndF = " F%.2f" % (60*speedExtrusion)
    cmndZ = " Z%.4f" % (Z)
    
    perimeters_command[i2] = perimeters_command[i2] + cmndG
    perimeters_command[i2] = perimeters_command[i2] + cmndF
    perimeters_command[i2] = perimeters_command[i2] + cmndX
    perimeters_command[i2] = perimeters_command[i2] + cmndY
    perimeters_command[i2] = perimeters_command[i2] + cmndZ
    perimeters_command[i2] = perimeters_command[i2] + "  ;Set Z"

    E = E + Ehop
    cmndG = "\nG1"
    cmndF = " F%.2f" % (60*speedExtrusion)
    cmndE = " E%.2f" % (E)
    cmndZ = " Z%.4f" % (Z)
    
    perimeters_command[i2] = perimeters_command[i2] + cmndG
    perimeters_command[i2] = perimeters_command[i2] + cmndF
    perimeters_command[i2] = perimeters_command[i2] + cmndE
    perimeters_command[i2] = perimeters_command[i2] + cmndZ
    perimeters_command[i2] = perimeters_command[i2] + "  ;Prime Nozzle"

    
    for i1 in list(range(1,len(airfoil[:,0]))):
        cmndG = "\nG1"
        
        nzlXo = perimeters[i1-1,0,i2]
        nzlYo = perimeters[i1-1,1,i2]
        nzlX = perimeters[i1,0,i2]
        nzlY = perimeters[i1,1,i2]
        nzlTrvl = np.sqrt((nzlX - nzlXo)**2 + (nzlY - nzlYo)**2)
        cmndX = " X%.4f" % (nzlX)
        cmndY = " Y%.4f" % (nzlY)
        
        if i2 == 0:
            Eprime = extrWidth*layer0Height
            dE = FrstLyrMult*FlowMult*Eprime*nzlTrvl/FilaDiameter**2
        else:
            Eprime = extrWidth*layerHeight
            dE = FlowMult*Eprime*nzlTrvl/FilaDiameter**2
            
        E = E + dE
        cmndE = " E%.4f" % (E)
        
        perimeters_command[i2] = perimeters_command[i2] + cmndG
        perimeters_command[i2] = perimeters_command[i2] + cmndF
        perimeters_command[i2] = perimeters_command[i2] + cmndE
        perimeters_command[i2] = perimeters_command[i2] + cmndX
        perimeters_command[i2] = perimeters_command[i2] + cmndY
        perimeters_command[i2] = perimeters_command[i2] + "  ;Perimeter"
    
    cmndG = "\nG1"
    
    nzlXo = perimeters[len(hull[:,0])-1,0,i2]
    nzlYo = perimeters[len(hull[:,0])-1,1,i2]
    nzlX = perimeters[0,0,i2]
    nzlY = perimeters[0,1,i2]
    
    nzlTrvl = np.sqrt((nzlX - nzlXo)**2 + (nzlY - nzlYo)**2)
    if i2 == 0:
        Eprime = extrWidth*layer0Height
        dE = FrstLyrMult*FlowMult*Eprime*nzlTrvl/FilaDiameter**2
    else:
        Eprime = extrWidth*layerHeight
        dE = FlowMult*Eprime*nzlTrvl/FilaDiameter**2
        
    E = E + dE
    cmndE = " E%.4f" % (E)
    
    cmndX = " X%.4f" % (perimeters[0,0,i2])
    cmndY = " Y%.4f" % (perimeters[0,1,i2])
    perimeters_command[i2] = perimeters_command[i2] + cmndG
    perimeters_command[i2] = perimeters_command[i2] + cmndF
    perimeters_command[i2] = perimeters_command[i2] + cmndE
    perimeters_command[i2] = perimeters_command[i2] + cmndX
    perimeters_command[i2] = perimeters_command[i2] + cmndY
    perimeters_command[i2] = perimeters_command[i2] + "  ;Perimeter End"
    
#Infill------------------------------------------------------------    
    cmndInfill = "\n"
    Infillx = c_airfoil*np.array(spar) + origin[0]
    Infilly = [[88,78.6],[87,80]]
    
    infillTrvl = [Infilly[0][0] - Infilly[0][1],Infilly[1][0] - Infilly[1][1]]
#     print(infillTrvl)
    
    E = E - Ehop
    cmndG = "G1"
    cmndE = " E%.4f" % (E)
    cmndZ = " Z%.4f" % (Z + Zhop)
    
    cmndInfill = cmndInfill + cmndG
    cmndInfill = cmndInfill + cmndE
    cmndInfill = cmndInfill + cmndZ
    cmndInfill = cmndInfill + "  ;Infill Hop"
    
    
    cmndG = "\nG0"
    cmndX = " X%.4f" % (Infillx[0])
    cmndY = " Y%.4f" % (Infilly[0][0])
    cmndZ = " Z%.4f" % (Z + Zhop)
    
    cmndInfill = cmndInfill + cmndG
    cmndInfill = cmndInfill + cmndX
    cmndInfill = cmndInfill + cmndY
    cmndInfill = cmndInfill + cmndZ
    cmndInfill = cmndInfill + "  ;Infill Position"
    
    cmndG = "\nG0"
    cmndX = " X%.4f" % (Infillx[0])
    cmndY = " Y%.4f" % (Infilly[0][0])
    cmndZ = " Z%.4f" % (Z)
    
    cmndInfill = cmndInfill + cmndG
    cmndInfill = cmndInfill + cmndX
    cmndInfill = cmndInfill + cmndY
    cmndInfill = cmndInfill + cmndZ
    cmndInfill = cmndInfill + "  ;Infill Position"
    
    E = E + Ehop
    cmndG = "\nG1"
    cmndE = " E%.4f" % (E)
    cmndZ = " Z%.4f" % (Z)
    
    cmndInfill = cmndInfill + cmndG
    cmndInfill = cmndInfill + cmndE
    cmndInfill = cmndInfill + cmndX
    cmndInfill = cmndInfill + cmndY
    cmndInfill = cmndInfill + cmndZ
    cmndInfill = cmndInfill + "  ;Infill Prime"
    
    cmndG = "\nG1"
    cmndX = " X%.4f" % (Infillx[0])
    cmndY = " Y%.4f" % (Infilly[0][1])
    
    if i2 == 0:
        Eprime = extrWidth*layer0Height
        dE = FrstLyrMult*FlowMult*Eprime*infillTrvl[0]/FilaDiameter**2
    else:
        Eprime = extrWidth*layerHeight
        dE = FlowMult*Eprime*infillTrvl[0]/FilaDiameter**2
    E = E + dE
    cmndE = " E%.4f" % (E)
    
    cmndInfill = cmndInfill + cmndG
    cmndInfill = cmndInfill + cmndF
    cmndInfill = cmndInfill + cmndE
    cmndInfill = cmndInfill + cmndX
    cmndInfill = cmndInfill + cmndY
    cmndInfill = cmndInfill + cmndZ
    cmndInfill = cmndInfill + "  ;Infill"
    
#----------------------------------------------------------------------    
    
    E = E - Ehop
    cmndG = "\nG1"
    cmndE = " E%.4f" % (E)
    cmndZ = " Z%.4f" % (Z + Zhop)
    
    cmndInfill = cmndInfill + cmndG
    cmndInfill = cmndInfill + cmndE
    cmndInfill = cmndInfill + cmndZ
    cmndInfill = cmndInfill + "  ;Infill Hop"
    
    
    cmndG = "\nG0"
    cmndX = " X%.4f" % (Infillx[1])
    cmndY = " Y%.4f" % (Infilly[1][0])
    cmndZ = " Z%.4f" % (Z + Zhop)
    
    cmndInfill = cmndInfill + cmndG
    cmndInfill = cmndInfill + cmndX
    cmndInfill = cmndInfill + cmndY
    cmndInfill = cmndInfill + cmndZ
    cmndInfill = cmndInfill + "  ;Infill Position"
    
    cmndG = "\nG0"
    cmndX = " X%.4f" % (Infillx[1])
    cmndY = " Y%.4f" % (Infilly[1][0])
    cmndZ = " Z%.4f" % (Z)
    
    cmndInfill = cmndInfill + cmndG
    cmndInfill = cmndInfill + cmndX
    cmndInfill = cmndInfill + cmndY
    cmndInfill = cmndInfill + cmndZ
    cmndInfill = cmndInfill + "  ;Infill Position"
    
    E = E + Ehop
    cmndG = "\nG1"
    cmndE = " E%.4f" % (E)
    cmndZ = " Z%.4f" % (Z)
    
    cmndInfill = cmndInfill + cmndG
    cmndInfill = cmndInfill + cmndE
    cmndInfill = cmndInfill + cmndX
    cmndInfill = cmndInfill + cmndY
    cmndInfill = cmndInfill + cmndZ
    cmndInfill = cmndInfill + "  ;Infill Prime"
    
    cmndG = "\nG1"
    cmndX = " X%.4f" % (Infillx[1])
    cmndY = " Y%.4f" % (Infilly[1][1])
    
    if i2 == 0:
        Eprime = extrWidth*layer0Height
        dE = FrstLyrMult*FlowMult*Eprime*infillTrvl[0]/FilaDiameter**2
    else:
        Eprime = extrWidth*layerHeight
        dE = FlowMult*Eprime*infillTrvl[0]/FilaDiameter**2
    E = E + dE
    cmndE = " E%.4f" % (E)
    
    cmndInfill = cmndInfill + cmndG
    cmndInfill = cmndInfill + cmndF
    cmndInfill = cmndInfill + cmndE
    cmndInfill = cmndInfill + cmndX
    cmndInfill = cmndInfill + cmndY
    cmndInfill = cmndInfill + cmndZ
    cmndInfill = cmndInfill + "  ;Infill"
    
    E = E - Ehop
    cmndG = "\nG1"
    cmndE = " E%.4f" % (E)
    cmndZ = " Z%.4f" % (Z + Zhop)
    
    cmndInfill = cmndInfill + cmndG
    cmndInfill = cmndInfill + cmndE
    cmndInfill = cmndInfill + cmndZ
    cmndInfill = cmndInfill + "  ;Infill Hop"
    
    perimeters_command[i2] = perimeters_command[i2] + cmndInfill
    
#----------------------------------------------------------------------
FilaLength = "Filament Used = %.3fmm" % (E)
PartMass = FilaDensity*(E/1000)*np.pi*(FilaDiameter/1000)**2/4
PartMassStr = "Filament Used = %.1fg" % (PartMass/1000)
PartCost = PartMass*FilaCost
PartCostStr = "Filament Cost = $%.2f" % (PartCost)


# ------------------------ Final GCode Write ------------------------- #
# Open a file
fo = open(fileGCodeOut+".gcode", "w")
fo.write(GCodeStart)
fo.write("\n\n")
fo.write(brim_command)
for i2 in list(range(n_lyrs)):
    fo.write(perimeters_command[i2])
fo.write("\n\n\n;"+FilaLength)
fo.write("\n;"+PartMassStr)
fo.write("\n;"+PartCostStr)
fo.write("\n\n\n")
fo.write(GCodeEnd)
fo.close() # Close opened file



print(FilaLength)
print(PartMassStr)
print(PartCostStr)
print("Time = ",time.time() - t)
fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.plot(hull[:,0],hull[:,1])
for i in list(range(skirtnum)):
    ax.plot(Skirts[:,0,i],Skirts[:,1,i])
ax.set_aspect(aspect=1)
plt.grid(True)
plt.show()

print("Done!")