def naca4digit(a,b,c,d):
	raise NotImplementedError
	
def plot_foil(foil):
    import numpy as np
    import matplotlib.pyplot as plt
    airfoil = np.array(foil)
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    ax.plot(airfoil[:,0],airfoil[:,1],marker='.')
    ax.set_aspect(aspect=1)
    plt.show()
    
def AFGeomAnalysis(AFcoords,nInterp):
    """
    Parameters
    ----------
    AFcoords : Selig-Style AF Coords
    nInterp : Scalar
        DESCRIPTION.

    Returns
    -------
    AFGeom : TYPE: numpy array
             DIMS: nInterp x 5
             Col0: x coords, range: [xLE,xTE], cosine discretization
             Col1: US y coords
             Col2: LS y coords
             Col3: thickness
             Col4: camber
    """
    import numpy as np
    from numpy import pi
    AFcoords = np.array(AFcoords)  # Convert into numpy array
    n = len(AFcoords)           # Measure number of points
    xTE = 0.5*(AFcoords[0,0] + AFcoords[-1,0])   # Find TE x-coord
    if AFcoords[n//2,0] < xTE: # xLE < xTE
        xLE = min(AFcoords[:,0])    # Find LE x-coord
    else: # Flying European Style!
        xLE = max(AFcoords[:,0])    # Find LE x-coord
    c = abs(xTE - xLE)
    
    # identify number of US points
    nUS = 0
    if xTE > xLE:
        for i in range(n-1):
            if AFcoords[i+1,0] < AFcoords[i,0]:
                nUS += 1
            else:
                break
    else:
        for i in range(n-1):
            if AFcoords[i+1,0] > AFcoords[i,0]:
                nUS += 1
            else:
                break
        
    nUS += 1
    USraw = np.zeros((nUS,2))
    nLS = n - nUS + 1
    LSraw = np.zeros((nLS,2))
    
    # create separate array for US points
    for i in range(nUS):
        USraw[i,:] = AFcoords[i,:]
    
    # create separate array for LS points
    for i in range(nLS):
        LSraw[i,:] = AFcoords[i+nUS-1,:]
    
    if xTE > xLE:
        USraw = np.flip(USraw,0)
    else:
        LSraw = np.flip(LSraw,0)
        
    AFGeom = np.zeros((nInterp,5)) # initialize output array
    dTheta = pi/(nInterp - 1)   # cosine discretization size
    
    # iterates from LE to TE to rediscretize airfoil
    for i in range(nInterp):
        AFGeom[i,0] = c*(0.5*np.cos(pi - i*dTheta) + 0.5) + min(xLE,xTE) #Col0 entry
        
        # interp ith US point
        for j in range(nUS - 1):
            if AFGeom[i,0] >= USraw[j,0] and AFGeom[i,0] < USraw[j+1,0]:
                AFGeom[i,1] = (((USraw[j+1,1] - USraw[j,1]) \
                /(USraw[j+1,0] - USraw[j,0])) \
                *(AFGeom[i,0] - USraw[j,0]) + USraw[j,1])
                break
        
        # interp ith LS point
        for j in range(nLS - 1):
            if AFGeom[i,0] >= LSraw[j,0] and AFGeom[i,0] < LSraw[j+1,0]:
                AFGeom[i,2] = (((LSraw[j+1,1] - LSraw[j,1]) \
                /(LSraw[j+1,0] - LSraw[j,0])) \
                *(AFGeom[i,0] - LSraw[j,0]) + LSraw[j,1])
                break
        
        AFGeom[i,3] =     (AFGeom[i,1] - AFGeom[i,2])   # ith thickness value
        AFGeom[i,4] = 0.5*(AFGeom[i,1] + AFGeom[i,2])   # ith camber value
    # for i in range(nInterp):
    #     AFGeom[i,0] = c*(0.5*np.cos(pi - i*dTheta) + 0.5) + xLE #Col0 entry
    # AFGeom[:,1] = np.interp(AFGeom[:,0],USraw[:,0],USraw[:,1] )
    # AFGeom[:,2] = np.interp(AFGeom[:,0],LSraw[:,0],LSraw[:,1] )
    # AFGeom[i,3] =     (AFGeom[i,1] - AFGeom[i,2])   # ith thickness value
    # AFGeom[i,4] = 0.5*(AFGeom[i,1] + AFGeom[i,2])   # ith camber value
    
    return AFGeom

def AFArc(AFcoords):
    """
    Parameters
    ----------
    AFcoords :  TYPE:        n x 2 numpy array-like
                DESCRIPTION: Selig-Style AF Coords

    Returns
    -------
    s : TYPE:        Scalar
        DESCRIPTION: Airfoil Perimeter
    """
    import numpy as np
    AFcoords = np.array(AFcoords)  # Convert into numpy array
    n = len(AFcoords)           # Measure number of points
    s = 0                           # Initialize perimeter at zero
    for i in range(n-1):
        dx = AFcoords[i+1,0]-AFcoords[i,0]
        dy = AFcoords[i+1,1]-AFcoords[i,1]
        ds = (dx**2 + dy**2)**0.5
        s = s + ds
    
    return s

def AFBendProperties(AFcoords,xMax, plotSwitch = 0):
    """
    Parameters
    ----------
    AFcoords :  TYPE:        n x 2 numpy array-like
                DESCRIPTION: Selig-Style AF Coords
    xMax :      TYPE:       Scalar
                Description: Max x/c location to consider AF skin for structure
                             (i.e. the flap/aileron hinge line)

    Returns
    -------
    IxxNorm :   TYPE:        Scalar
                DESCRIPTION: Ixx/(t_skin*c^3)
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from CodeFiles.Aerodynamics.Airfoils.operations import AFGeomAnalysis
    AFcoords = np.array(AFcoords)        # Convert into numpy array
    
    # Convert AF to [x,yUS,yLS,t,h]
    nInterp  = 1001
    AFGeom   = AFGeomAnalysis(AFcoords,nInterp)
    AFGeomTrimmed  = AFGeom[AFGeom[:,0]<=xMax,:] # only consider foreward of hinge
    
    # breakpoint()
    
    nAFGeomTrimmed = len(AFGeomTrimmed)
    
    xTau = np.argmax(AFGeomTrimmed[:,3])/nInterp - 0.07
    dx = 30/122
    xSpars  = [xTau - dx, xTau, xTau + dx, xTau + 2*dx]
    USspars = np.interp(xSpars, AFGeom[:,0], AFGeom[:,1])
    LSspars = np.interp(xSpars, AFGeom[:,0], AFGeom[:,2])
    
    # breakpoint()
    
    # calculate  centroid [xbar,ybar]
    sum_xds = 0                          # Initialize
    sum_yds = 0                          # Initialize
    sum_ds  = 0                          # Initialize
    for i in range(nAFGeomTrimmed - 1): #iterate from LE to xHinge
        x    = 0.5*(AFGeomTrimmed[i+1,0] + AFGeomTrimmed[i,0])
        dx   = AFGeomTrimmed[i+1,0] - AFGeomTrimmed[i,0]
        yUS  = 0.5*(AFGeomTrimmed[i+1,1] + AFGeomTrimmed[i,1])
        yLS  = 0.5*(AFGeomTrimmed[i+1,2] + AFGeomTrimmed[i,2])
        dyUS = AFGeomTrimmed[i+1,1] - AFGeomTrimmed[i,1]
        dyLS = AFGeomTrimmed[i+1,2] - AFGeomTrimmed[i,2]
        dsUS = (dx**2 + dyUS**2)**0.5
        dsLS = (dx**2 + dyLS**2)**0.5
        
        sum_ds  = sum_ds  +     dsUS +     dsLS
        sum_xds = sum_xds +   x*dsUS +   x*dsLS
        sum_yds = sum_yds + yUS*dsUS + yLS*dsLS
    # Handle blunt hinge
    x   = AFGeomTrimmed[-1,0]
    yUS = AFGeomTrimmed[-1,1]
    yLS = AFGeomTrimmed[-1,2]
    y   = 0.5*(yUS + yLS)
    ds  = yUS - yLS  # no dx, only consider dy
    sum_ds  = sum_ds  +   ds
    sum_xds = sum_xds + x*ds
    sum_yds = sum_yds + y*ds
    
    # calculate xbar and ybar
    xbar = sum_xds/sum_ds
    ybar = sum_yds/sum_ds
    
    # calculate bending intertia [IxxNorm, IyyNorm]
    sum_Ixx = 0
    sum_Iyy = 0
    for i in range(nAFGeomTrimmed - 1):
        x    = 0.5*(AFGeomTrimmed[i+1,0] + AFGeomTrimmed[i,0])
        yUS  = 0.5*(AFGeomTrimmed[i+1,1] + AFGeomTrimmed[i,1])
        yLS  = 0.5*(AFGeomTrimmed[i+1,2] + AFGeomTrimmed[i,2])
        dx   = AFGeomTrimmed[i+1,0] - AFGeomTrimmed[i,0]
        dyUS = AFGeomTrimmed[i+1,1] - AFGeomTrimmed[i,1]
        dyLS = AFGeomTrimmed[i+1,2] - AFGeomTrimmed[i,2]
        dsUS = (dx**2 + dyUS**2)**0.5
        dsLS = (dx**2 + dyLS**2)**0.5
        
        sum_Ixx = sum_Ixx + (yUS-ybar)**2*dsUS + (yLS-ybar)**2*dsLS
        sum_Iyy = sum_Iyy + (x-xbar)**2*(dsUS + dsLS)
    
    tau_Hinge = AFGeomTrimmed[-1,3]
    eta_Hinge = AFGeomTrimmed[-1,4]
    IxxNorm = sum_Ixx + (1/12)*tau_Hinge**3 + tau_Hinge*(eta_Hinge - ybar)**2
    IyyNorm = sum_Iyy + (1/12)*tau_Hinge +    tau_Hinge*(xMax - xbar)**2
    JzNorm = IxxNorm + IyyNorm
    
    if plotSwitch:
        width  = 6.5
        height = 7
        tauMax = max(AFGeom[:,3])
        ylims = [ybar - 2*tauMax, ybar + 2*tauMax] 
        plt.figure(figsize=(width, height))
        # plt.subplot(211)
        plt.grid()
        plt.plot(AFGeom[:,0],AFGeom[:,1],'b:')
        plt.plot(AFGeom[:,0],AFGeom[:,2],'r:')
        # plt.plot(AFGeom[:,0],AFGeom[:,4],'k:')
        plt.plot([-1,2],[ybar,ybar],'k:')
        plt.plot([xbar,xbar],[-1,1],'k:')
        plt.plot(AFGeomTrimmed[:,0],AFGeomTrimmed[:,1],'b-')
        plt.plot([xMax,xMax],[AFGeomTrimmed[-1,1],eta_Hinge],'b-')
        plt.plot(AFGeomTrimmed[:,0],AFGeomTrimmed[:,2],'r-')
        plt.plot([xMax,xMax],[AFGeomTrimmed[-1,2],eta_Hinge],'r-')
        # for i in range(nSpars):
        #     xCoords = [sparMat[i,0], sparMat[i,0]]
        #     yCoords = [sparMat[i,1], sparMat[i,2]]
        #     plt.plot(xCoords, yCoords,'k-')
        plt.gca().set_aspect("equal")
        plt.xlim( -0.05   , 1.05 )
        plt.ylim(ylims)
        plt.xticks(fontsize = 16)
        plt.yticks(fontsize = 16)
        plt.ylabel(r'$y/c$', fontsize=20)
        plt.xlabel(r'$x/c$', fontsize=20)
        titleStr = """$\overline{x}/c$ = %4.3f\n
$\overline{y}/c$ = %4.3f\n
$I_{x}/(t_{sk}*c^3)$ = %5.4e\n
$J_{z}/(t_{sk}*c^3)$ = %5.4e""" %(xbar, ybar, IxxNorm, JzNorm)
        plt.title(titleStr, fontsize=20)
        
        
        c = 122
        width  = 6.5
        height = 3
        tauMax = max(AFGeom[:,3])
        ylims = [c*(ybar - 2*tauMax), c*(ybar + 2*tauMax)] 
        plt.figure(figsize=(width, height))
        # plt.subplot(211)
        plt.grid()
        plt.plot(c*AFGeom[:,0],c*AFGeom[:,1],'k:')
        plt.plot(c*AFGeom[:,0],c*AFGeom[:,2],'k:')
        # plt.plot(AFGeom[:,0],AFGeom[:,4],'k:')
        # plt.plot([-1,2],[ybar,ybar],'k:')
        # plt.plot([xbar,xbar],[-1,1],'k:')
        plt.plot(c*AFGeomTrimmed[:,0],c*AFGeomTrimmed[:,1],'k-')
        plt.plot([c*xMax,c*xMax],[c*AFGeomTrimmed[-1,1],c*eta_Hinge],'k-')
        plt.plot(c*AFGeomTrimmed[:,0],c*AFGeomTrimmed[:,2],'k-')
        plt.plot([c*xMax,c*xMax],[c*AFGeomTrimmed[-1,2],c*eta_Hinge],'k-')
        for i in range(4):
            xCoords = [c* xSpars[i],c* xSpars[i]]
            yCoords = [c*USspars[i],c*LSspars[i]]
            plt.plot(xCoords, yCoords,'k-')
        plt.gca().set_aspect("equal")
        plt.xlim( -c*0.05   , c*1.05 )
        plt.ylim(ylims)
        plt.xticks(fontsize = 16)
        plt.yticks(fontsize = 16)
        plt.ylabel(r'$y, mm$', fontsize=20)
        plt.xlabel(r'$x, mm$', fontsize=20)
        # plt.title(titleStr, fontsize=20)

    return [IxxNorm,JzNorm,ybar]

def AFCurveAnalysis(AFcoords,nInterp):
    """
    Parameters
    ----------
    AFcoords : Selig-Style AF Coords
    nInterp : Scalar
    Description: Calculates curvature on the upper and lower surface of the airfoil

    Returns
    -------
    AFcurve : TYPE: numpy array
             DIMS: nInterp x 3
             Col0: x coords, range: [xLE,xTE], cosine discretization
             Col1: US curvature
             Col2: LS curvature
    """
    import numpy as np
    from numpy import pi
    import copy
    AFcoords = np.array(AFcoords)  # Convert into numpy array
    n = len(AFcoords)           # Measure number of points
    xTE = 0.5*(AFcoords[0,0] + AFcoords[-1,0])   # Find TE x-coord
    if AFcoords[n//2,0] < xTE: # xLE < xTE
        xLE = min(AFcoords[:,0])    # Find LE x-coord
    else: # Flying European Style!
        xLE = max(AFcoords[:,0])    # Find LE x-coord
    c = abs(xTE - xLE)
    
    # identify number of US points
    nUS = 0
    if xTE > xLE:
        for i in range(n-1):
            if AFcoords[i+1,0] < AFcoords[i,0]:
                nUS += 1
            else:
                break
    else:
        for i in range(n-1):
            if AFcoords[i+1,0] > AFcoords[i,0]:
                nUS += 1
            else:
                break
        
    nUS += 1
    USraw = np.zeros((nUS,2))
    nLS = n - nUS + 1
    LSraw = np.zeros((nLS,2))
    
    # create separate array for US points
    for i in range(nUS):
        USraw[i,:] = AFcoords[i,:]
    
    # create separate array for LS points
    for i in range(nLS):
        LSraw[i,:] = AFcoords[i+nUS-1,:]
    
    if xTE > xLE:
        USraw = np.flip(USraw,0)
    else:
        LSraw = np.flip(LSraw,0)
    
    # Take required gradients of airfoil
    dydxUS = np.gradient(USraw[:,1],USraw[:,0],edge_order=2)
    dydxLS = np.gradient(LSraw[:,1],LSraw[:,0],edge_order=2)
    d2ydx2US = np.gradient(dydxUS,USraw[:,0],edge_order=2)
    d2ydx2LS = np.gradient(dydxLS,LSraw[:,0],edge_order=2)
    
    # calculate curvature
    UScurve =  copy.deepcopy(USraw)
    LScurve =  copy.deepcopy(LSraw)
    for i in range(len(USraw[:,0])):
        UScurve[i,1] = abs(d2ydx2US[i])/((1+dydxUS[i]**2)**(3/2))
    for i in range(len(LSraw[:,0])):
        LScurve[i,1] = abs(d2ydx2LS[i])/((1+dydxLS[i]**2)**(3/2))
    
    dTheta = pi/(nInterp - 1)   # cosine discretization size
    AFcurve = np.zeros((nInterp,3))  # initialize output array
    for i in range(nInterp): # determine cosine spacing
        AFcurve[i,0] = c*(0.5*np.cos(pi - i*dTheta) + 0.5) + min(xLE,xTE) #Col0 entry
    
    # Interpolate US and LS curvature to desired spacing
    AFcurve[:,1] = np.interp(AFcurve[:,0],UScurve[:,0],UScurve[:,1])
    AFcurve[:,2] = np.interp(AFcurve[:,0],LScurve[:,0],LScurve[:,1])

    return AFcurve

def CircProfile(k_US,k_LS,chord,n):
    import numpy as np
    if k_US == 0:
        k_US = 0.000001
    if k_LS == 0:
        k_LS = 0.000001

    rUS = 1/k_US
    rLS = 1/k_LS
    rLE = 0.05*chord
    
    if 2*rUS < chord: 
        print("Upper Surface too curved!!!")
        return(0)
    if 2*rLS < chord: 
        print("Lower Surface too curved!!!")
        return(0)

    nUS = int(0.42*n)
    nLS = int(0.42*n)
    nLE = n - nUS - nLS

    circUS = np.zeros((nUS,2))
    circLS = np.zeros((nLS,2))
    circLE = np.zeros((nLE,2))

    # Generate upper surface
    thetaUS  = np.arccos(((rUS - rLE)**2 + rUS**2 - (chord - rLE)**2)/(2*rUS*(rUS-rLE)))
    phiUS    = np.arccos(((chord - rLE)**2 + (rUS - rLE)**2 - rUS**2)/(2*(chord - rLE)*(rUS - rLE)))
    thetaUSo = np.pi - phiUS - thetaUS
    cUS      = [rLE + (rUS - rLE)*np.cos(phiUS), -(rUS - rLE)*np.sin(phiUS)]

    dTheta = thetaUS/nUS
    for i in range(nUS):
        xnew = cUS[0] + rUS*np.cos(i*dTheta + thetaUSo)
        ynew = cUS[1] + rUS*np.sin(i*dTheta + thetaUSo)
        circUS[i,:] = [xnew,ynew]

    # Generate lower surface
    thetaLS  = np.arccos(((rLS + rLE)**2 + rLS**2 - (chord - rLE)**2)/(2*rLS*(rLS + rLE)))
    phiLS    = np.arccos(((chord - rLE)**2 + (rLS + rLE)**2 - rLS**2)/(2*(chord - rLE)*(rLS + rLE)))
    thetaLSo = np.pi - phiLS - thetaLS
    cLS      = [rLE + (rLS + rLE)*np.cos(phiLS), -(rLS + rLE)*np.sin(phiLS)]
    # print(phiUS*180/np.pi)

    dTheta = thetaLS/nLS
    for i in range(nLS):
        xnew = cLS[0] + rLS*np.cos(i*dTheta + thetaLSo)
        ynew = cLS[1] + rLS*np.sin(i*dTheta + thetaLSo)
        circLS[i,:] = [xnew,ynew]
    circLS = np.flip(circLS,axis = 0)

    gamma0 = np.pi - phiUS
    gamma1 = 2*np.pi - gamma0 - phiLS
    dGamma = gamma1/(nLE - 1)
    for i in range(nLE):
        circLE[i,0] = rLE*np.cos(i*dGamma + gamma0)
        circLE[i,1] = rLE*np.sin(i*dGamma + gamma0)

    circLE = circLE + np.array([rLE,0])

    profile = np.concatenate((circUS, circLE, circLS))
    profile = profile/chord
    theta = np.arcsin(rLE/(chord - rLE))
    c = np.cos(theta)
    s = np.sin(theta)
    rotMat = np.matrix([[ c,-s],[ s, c]])
    profile = profile@rotMat
    profile[:,1] = profile[:,1] - min(profile[:,1])
    
    return(profile)

def plotAFCurve(chord,AFcoords,SparLoc):
    import matplotlib.pyplot as plt
    import numpy as np
    from CodeFiles.Aerodynamics.Airfoils.operations import AFGeomAnalysis
    from CodeFiles.Aerodynamics.Airfoils.operations import AFCurveAnalysis
    AFcoordsMat = np.array(AFcoords)
    
    AFGeom  = AFGeomAnalysis(AFcoordsMat,100)
    AFCurve = AFCurveAnalysis(chord*AFcoordsMat,100)
    
    USCurveSample = AFCurve[50,1]
    
    nSpars = len(SparLoc)
    sparMat = np.zeros((nSpars,3))
    sparMat[:,0] = SparLoc
    sparMat[:,1] = np.interp(SparLoc, AFGeom[:,0], AFGeom[:,1])
    sparMat[:,2] = np.interp(SparLoc, AFGeom[:,0], AFGeom[:,2])
    dx = round(1000*chord*(SparLoc[1]-SparLoc[0]))
    
    width  = 3.25
    height = 3.5
    
    plt.figure(figsize=(width, height))
    plt.subplot(211)
    plt.grid()
    plt.plot(AFGeom[:,0],AFGeom[:,1],'b-')
    plt.plot(AFGeom[:,0],AFGeom[:,2],'r-')
    for i in range(nSpars):
        xCoords = [sparMat[i,0], sparMat[i,0]]
        yCoords = [sparMat[i,1], sparMat[i,2]]
        plt.plot(xCoords, yCoords,'k-')
    plt.gca().set_aspect("equal")
    plt.xlim( 0   , 1 )
    plt.ylim(-0.20, 0.25 )
    plt.ylabel(r'$y/c$', fontsize=10)
    plt.title('Chord = '+str(1000*chord)+"mm; $dx$ = "+str(dx)+"mm", fontsize=10)
    
    plt.subplot(212)
    plt.grid()
    plt.plot(AFGeom[:,0],AFCurve[:,1],'b-')
    plt.plot(AFGeom[:,0],AFCurve[:,2],'r-')
    # plt.gca().set_aspect("equal")
    plt.xlim( 0   , 1)
    plt.ylim(0, 2*round(USCurveSample) )
    plt.xlabel(r'$x/c$', fontsize=10)
    plt.ylabel(r'$\kappa$[m]', fontsize=10)
    
    
    width  = 2
    height = 1.5
    
    plt.figure(figsize=(width, height))
    plt.grid()
    plt.plot(AFGeom[:,0],AFGeom[:,1],'k-')
    plt.plot(AFGeom[:,0],AFGeom[:,2],'k-')
    for i in range(nSpars):
        xCoords = [sparMat[i,0], sparMat[i,0]]
        yCoords = [sparMat[i,1], sparMat[i,2]]
        plt.plot(xCoords, yCoords,'k-')
    plt.gca().set_aspect("equal")
    plt.xlim( -0.1   , 1.1 )
    plt.ylim(-0.10, 0.25 )
    plt.ylabel(r'$y/c$', fontsize=10)
    plt.xlabel(r'$x/c$', fontsize=10)
    # plt.title('Chord = '+str(1000*chord)+"mm; $dx$ = "+str(dx)+"mm", fontsize=10)
    