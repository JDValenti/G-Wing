# -*- coding: utf-8 -*-

def FuseWeightDistCalc(P,b,n,fuseY,fuseW,y_vec):
    import numpy as np
    pFuse_vec = np.zeros((n,1))
    fuseCount = 0
    # fusei = 
    for i in range(n):
        if y_vec[i] >= (fuseY - 0.5*fuseW) and y_vec[i] <= (fuseY + 0.5*fuseW):
            fuseCount += 1
            fuseEnd = i
    if fuseCount < 2:
        fusei = n//2
        if n%2 == 1: # Define p if n is odd
            pFuse_vec[fusei] = (n/b)*P
        else:        # Define p if n is even
            pFuse_vec[fusei - 1] = 0.5*(n/b)*P
            pFuse_vec[fusei]     = 0.5*(n/b)*P
    else:
        fuseStart = fuseEnd - fuseCount
        pFuse_vec[fuseStart+1:fuseEnd+1] = (1/fuseCount)*(n/b)*P*np.ones((fuseCount,1))
        fusei = fuseEnd - fuseCount//2
    return(pFuse_vec)

def Forces_n_Moments(b,yb_vec,l_vec,p_vecs,wStruct_vec,dispPlots,tipLoaded):
    import numpy as np   
    n = len(l_vec)
    e_vec = np.ones((n,1))
    # Defining K matrix
    K = np.zeros((n,n))        # initialize K matrix
    for k in range(n):      # nest for loops assign elements
        for j in range(n):
            # K[j,k] = 0.5*b*(1/n)**2*((j-k)**2)**0.5
            K[j,k] = (1/2)*(b**2/n**2)*np.abs(j-k)
            # K[j,k] = b*(1/n)**2*((j-k)**2)**0.5
    import copy
    l_vec_Temp = copy.deepcopy(l_vec)
    l_vec = np.zeros((n,1))
    for j in range(n):
        l_vec[j,0] = l_vec_Temp[j]
    del l_vec_Temp
    if tipLoaded:
        tipLoad = sum(l_vec)[0]/2
        l_vec = np.concatenate([[[tipLoad]],np.zeros((n-2,1)),[[tipLoad]]])
    n_pvecs = len(p_vecs[0,:])
    wTot_vec = p_vecs@np.ones((n_pvecs,1)) + wStruct_vec # Summing weights
    f_vec = l_vec - wTot_vec   # Summing forces
    mb_vec = K@f_vec   # Calculating Bending Moments
    
    # Analytical Comparison
    L_tot = ((b/n)*e_vec.T@l_vec)[0,0]
    W_tot = ((b/n)*e_vec.T@wTot_vec)[0,0]
    f_net = (e_vec.T@f_vec)[0,0]
    print("Total Lift   = ",1000*L_tot/9.8,"g")
    print("Total Weight = ",1000*W_tot/9.8,"g")
    print("Net Force    = ",1000*f_net/9.8,"g")
    
    bSemi = b/2
    m_root_rect = (0.5*L_tot)*(0.50*bSemi)  # Assume Constant Lift Dist
    m_root_tri  = (0.5*L_tot)*(0.33*bSemi)  # Assume Tianglular Lift Dist
    
    l_const_vec = (L_tot/b)*np.ones((n,1))
    m_const_vec = K@(l_const_vec - wTot_vec)
    
    if dispPlots:
        import matplotlib.pyplot as plt
        # width  = 3.25
        # height = 3.5
        width  = 6.5
        height = 6.6
        plt.figure(figsize=(width, height))
        
        plt.subplot(211)
        plt.plot(yb_vec,l_vec,'b--',label=r'Lift',linewidth=2)
        for i in range(n_pvecs):
            if i == 0:
                lbl = 'Fuselage \nWeight'
                style = 'r-.'
            else:
                lbl = 'Wing \nWeight'
                style = 'r:'
            plt.plot(yb_vec,p_vecs[:,i], style, label=lbl, linewidth=2)
        # plt.plot(yb_vec,f_vec,'k-',label=r'Net',linewidth=2)
        # plt.plot(yb_vec,l_const_vec,'r',label=r'const',linewidth=2)
        # plt.plot(yb_vec,p_vecs[:,0], linestyle=':',label=r'Fuselage Weight',linewidth=1)
        # plt.plot(yb_vec,p_vecs[:,1], linestyle=':',label=r'Wing Weight',linewidth=1)
        # plt.xlabel("y");
        plt.ylabel("Force per\n Unit Span, N/m", fontsize=20);
        plt.yticks(fontsize=18);
        plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)), ());
        plt.legend(fontsize=18,loc='lower center', bbox_to_anchor=(0.5, 1.01),ncol = 3);
        plt.xlim(-0.5,0.5)
        # plt.ylim(top = 3*max(p_vecs[:,0]))
        plt.grid()    
        
        plt.subplot(212)
        plt.plot(yb_vec,mb_vec,label=r'Bending Moment', color='k', linewidth=2)
        # plt.plot(yb_vec,m_const_vec,label=r'Bending Moment', color='r', linewidth=2)
        # plt.plot(0,m_root_rect, 'rv')
        # plt.plot(0,m_root_tri, 'b^')
        plt.xlabel("y");
        plt.ylabel('Bending \n'+r'Moment, N $\times$ m', fontsize=20);
        plt.yticks(fontsize=18);
        plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)), fontsize=18);
        plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
        plt.xlim(-0.5,0.5)
        plt.grid()
        
    return(mb_vec)

def sectionSparPlace(m,mMax,xtmax,dxLim):
    """
    This function generates the chordwise location of all spars for a given 
    section based on a scalars m and dxlim
    """
    import numpy as np
    mMax = int(np.ceil(mMax))
    # m = int(m)

    # Method 1 **********************************
    # Initialize x/c values of spars [Main, Fore1, Aft1, Fore2, Aft2... ]
    try:
        xcSpars  = np.zeros((1 , 2*mMax + 1))
    except:
        breakpoint()
    # breakpoint()
    xcSpars[0,  0] = xtmax
    
    dxFor =    xtmax /(m+1)
    dxAft = (1-xtmax)/(m+1)
    for i in range(mMax):
        # breakpoint()
        xcSpars[0,  2*i + 1] = max(xtmax - (i+1)*dxFor,0)
        xcSpars[0,2*(i + 1)] = min(xtmax + (i+1)*dxAft,1)

    # #  Method 2 *****************************
    # xcSpars = [xtmax]
    # dxFor = min(   xtmax /(m+1),dxLim)
    # for i in range(50):
    #     xc = xtmax - (i+1)*dxFor
    #     if xc > 0:
    #         xcSpars.append(xc)
    #     else:
    #         break
        
    # dxAft = min((1-xtmax)/(m+1),dxLim)
    
    # for i in range(50):
    #     xc = xtmax + (i+1)*dxAft
    #     if xc < 1:
    #         xcSpars.append(xc)
    #     else:
    #         break
    # for i in range(len(xcSpars)):
    #     xcSpars[i] = float(xcSpars[i])
    
    # xcSpars.sort()
    
    return(xcSpars)

# def wingSparNormPlace(m_vec,xtmax,dxLim_vec):
#     import numpy as np
#     n = len(m_vec)
#     xcSpars = ['']*n
#     nSpars = ['']*n
#     for j in range(n):
#         xcSpars[j] = sectionSparPlace(m_vec[j][0],int(np.ceil(max(m_vec))),xtmax,dxLim_vec[j])
#         nSpars[j]  = len(xcSpars[j])
    
#     breakpoint()
#     return(xcSpars)

def wingSparAbsPlace(xcSpars,c_vec,xcLL):
    import numpy as np
    n = len(c_vec)
    try:
        nSpars = len(xcSpars[0,:])
        xSpars = np.zeros((n,nSpars))
        for j in range(n):
            xSpars[j,:]= (xcSpars[j,:] - xcLL)*c_vec[j]
    except:
        nSpars = 1
        xSpars = np.zeros((n,nSpars))
        xSpars= (xcSpars - xcLL)*c_vec
    
    return(xSpars)

def wingSparAbs2NormPlace(xSpars,c_vec,xcLL):
    import numpy as np
    n = len(c_vec)
    if n != 1:
        xcSpars = np.zeros(np.shape(xSpars))
        for j in range(n):
            xcSpars[j,:]= xSpars[j,:]/c_vec[j] + xcLL
    else:
        xcSpars= xSpars/c_vec + xcLL
    return(xcSpars)

def wingSparPlace(m_vec,xtmax,c_vec,xcLL,dxLim):
    
    import numpy as np
    n = len(m_vec)
    # xcSpars = ['']*n
    # xSpars  = ['']*n
    # for j in range(n):     
    #     # breakpoint()
    #     xcSpars[j] = sectionSparPlace(m_vec[j],xtmax,dxLim/c_vec[j,0])
        
    #     nSpars = len(xcSpars[j])
    #     xSpars[j] = ['']*nSpars
    #     for k in range(nSpars):
    #         xSpars[j][k]  = (xcSpars[j][k] - xcLL)*c_vec[j][0]
    #     # breakpoint()
    mCeilMax = int(np.ceil(max(m_vec)))
    nSparsMax = 2*mCeilMax + 1
    xcSpars = np.zeros((n,nSparsMax))
    xSpars  = np.zeros((n,nSparsMax))
    for i in range(mCeilMax+1):
        xcSpars[:,2*i] = np.ones(n)
    for j in range(n):
        # nSpars = 2*int(np.ceil(m_vec[j]))+1
        
        xcSpars[j,:] = sectionSparPlace(m_vec[j],mCeilMax,xtmax,dxLim)
        
        for k in range(nSparsMax):
            xSpars[j,k]  = (xcSpars[j,k] - xcLL)*c_vec[j][0]
        

    arc     = np.zeros((n,2))
    for j in range(n):
        arc[j,0]     =          min(c_vec[j]*xtmax/(m_vec[j] + 1),dxLim)
        arc[j,1]     = min(c_vec[j]*min((1 - xtmax)/(m_vec[j] + 1),(xcLL - xtmax)),dxLim)
    
    # breakpoint()
    return(xcSpars,xSpars,arc)
    

def rSpars(c_vec,xcSpars,tauDist):
    """
    This function returns a spanwise distibution vector 
    that represents the sum of the the spar heights at 
    any given spanwise location.  
    This can easily be used to generate wing weight: 
        ww_vec = rho_m*g*w_e*rSpars_vec 
    The rSpars_vec elements can be summed to make the
        make the cost function to be minimized
    """
    import numpy as np
    n = len(c_vec)
    # nSpars = len(xcSpars[0,:])
    rSpars_vec = np.zeros((n,1))
    for j in range(n):
        # breakpoint()
        tauSpars = np.interp(xcSpars[j],tauDist[:,0],tauDist[:,1])
        rSpars_vec[j,0] = c_vec[j,0]*sum(tauSpars)
    return(rSpars_vec)

def bendingInertia(c_vec,xcSpars,tauDist,w_e):
    """
    Calculates the bending intertias of the spars 
    and returns a numpy array of the same shape
    as xcSpars containing the 
    """
    import numpy as np
    n = len(c_vec)
    import copy
    IxSpars = copy.deepcopy(xcSpars)
    for j in range(n):
        tauSpars = np.interp(xcSpars[j,:],tauDist[:,0],tauDist[:,1])
        for i in range(len(tauSpars)):
            IxSpars[j,i] = (1/12)*w_e*(tauSpars[i]*c_vec[j,0])**3
    return(IxSpars)

def bendingInertiaWAil(c_vec,xcSpars,tauDist,w_e,xc_hinge):
    """
    Calculates the bending intertias of the spars 
    and returns a numpy array of the same shape
    as xcSpars containing the 
    """
    import numpy as np
    # breakpoint()
    n = len(c_vec)
    import copy
    IxSpars = copy.deepcopy(xcSpars)
    for j in range(n):
        tauSpars = np.interp(xcSpars[j],tauDist[:,0],tauDist[:,1])
        for i in range(len(tauSpars)):
            try:
                if xcSpars[j,i] < xc_hinge:
                    IxSpars[j,i] = (1/12)*w_e*(tauSpars[i]*c_vec[j,0])**3
                else:
                    IxSpars[j,i] = 0
            except:
                breakpoint()
    # breakpoint()
    return(IxSpars)

def sectionBendingInertia(c,xcSpars,tauDist,w_e):
    """
    Calculates the bending intertias of the spars 
    and returns a numpy array of the same shape
    as xcSpars containing the 
    """
    import numpy as np
    IxSpars = 0
    nSpars  = np.size(xcSpars)
    tauSpars = np.interp(xcSpars,tauDist[:,0],tauDist[:,1])
    for i in range(nSpars):
        IxSpars = (1/12)*w_e*(tauSpars[0,i]*c)**3 + IxSpars
    return(IxSpars)

def sectionBendingInertiaWAil(c,xcSpars,tauDist,w_e,xc_hinge):
    """
    Calculates the bending intertias of the spars 
    and returns a numpy array of the same shape
    as xcSpars containing the 
    """
    import numpy as np
    IxSpars = 0
    nSpars  = np.size(xcSpars)
    tauSpars = np.interp(xcSpars,tauDist[:,0],tauDist[:,1])
    # breakpoint()
    # print("nSpars  = ", nSpars)
    # print("xcSpars = ", xcSpars)
    for i in range(nSpars):
        # print(i)
        # if i == 1:
        #     breakpoint()
        if xcSpars[0][i] < xc_hinge:
            
            IxSpars = (1/12)*w_e*(tauSpars[0][i]*c)**3 + IxSpars
    return(IxSpars)

def bisection(f,a,b,N):
    '''
    Function adapted from:
    https://www.math.ubc.ca/~pwalls/math-python/roots-optimization/bisection/
    
    Approximate solution of f(x)=0 on interval [a,b] by bisection method.
    Parameters
    ----------
    f : function
        The function for which we are trying to approximate a solution f(x)=0.
    a,b : numbers
        The interval in which to search for a solution. The function returns
        None if f(a)*f(b) >= 0 since a solution is not guaranteed.
    N : (positive) integer
        The number of iterations to implement.

    Returns
    -------
    x_N : number
        The midpoint of the Nth interval computed by the bisection method. The
        initial interval [a_0,b_0] is given by [a,b]. If f(m_n) == 0 for some
        midpoint m_n = (a_n + b_n)/2, then the function returns this solution.
        If all signs of values f(a_n), f(b_n) and f(m_n) are the same at any
        iteration, the bisection method fails and return None.

    Examples
    --------
    >>> f = lambda x: x**2 - x - 1
    >>> bisection(f,1,2,25)
    1.618033990263939
    >>> f = lambda x: (2*x - 1)*(x - 3)
    >>> bisection(f,0,1,10)
    0.5
    '''
    if f(a) >= 0:
        # print("MAIN SPAR STRONG ENOUGH")
        return 0
    if f(a)*f(b) >= 0:
        print("BISECTION METHOD FAILED!!! f(a) & f(b) ARE SAME SIGN!!!")
        return None
    a_n = a
    b_n = b
    for n in range(1,N+1):
        m_n = (a_n + b_n)/2
        f_m_n = f(m_n)
        if f(a_n)*f_m_n < 0:
            a_n = a_n
            b_n = m_n
        elif f(b_n)*f_m_n < 0:
            a_n = m_n
            b_n = b_n
        elif f_m_n == 0:
            print("Found exact solution.")
            return m_n
        else:
            print("BISECTION METHOD FAILED!!! BOUNDS DON'T STRADDLE ZERO")
            return None
    return (a_n + b_n)/2

def SkinWeightCalc(c_vec,AFcoords,b,w_e,rho_m,g):
    """
    Calculates the skin weight
    """
    import numpy as np
    from CodeFiles.Aerodynamics.Airfoils.operations import AFArc
    # b = 0.001*b
    # w_e = 0.001*w_e
    sc = AFArc(AFcoords) # Perimeter/chord of airfoil
    s_vec = sc*b*c_vec   # AF perimeters across the span
    n = len(c_vec)
    w_sk_array = rho_m*g*w_e*s_vec
    w_sk_vec = np.zeros((n,1))
    for j in range(n):
        w_sk_vec[j,0] = w_sk_array[j]
    return(w_sk_vec)

def RibLocations(b,yMax,ailStart):
    """
    Returns array of normalized rib locations (y/b)
    """
    import numpy as np
    nSemiSections = int(np.ceil((b/2 - ailStart)/yMax))
    sectSpan = (b/2 - ailStart)/nSemiSections
    ribLocSW = np.zeros(nSemiSections)
    for i in range(nSemiSections):
        ribLocSW[i] = ailStart + i*sectSpan
    ribLoc = np.concatenate((np.flip(-ribLocSW),ribLocSW))/b
    return(ribLoc)

def WingBend(b,E,I_vec,l_vec,p_vec,ww_vec,dispPlots):
    import numpy as np
    import matplotlib.pyplot as plt
    
    try:
        n =len(I_vec)
        
        dy = b/n
        L = sum(l_vec)*dy
        
        print("Lift for Bending          = %8.4f N = %8.2f g" %(L, 1000*L/9.8))
        wLoad = 0.5*L
        print("Sandbag Weight from Lift  = %8.4f N = %8.2f g" %(wLoad, 1000*wLoad/9.8))
        
        if 1: 
            wLoad = 243*(9.8/1000)
            print("User Input Sandbag Weight = %8.4f N = %8.2f g" %(wLoad, 1000*wLoad/9.8))
        
        
        y = np.linspace(-0.5*b,0.5*b,n)
        
        # y_canti = y[p_vec[:,0]==0]
        # y_canti = y_canti[len(y_canti)//2:len(y_canti)]
        # b_canti = max(y_canti) - min(y_canti)
        # l_canti = l_vec[p_vec[:,0]==0]*dy
        # l_canti = l_canti[len(l_canti)//2:len(l_canti)]
        # ww_canti = (p_vec[:,1]+ww_vec[:,0])[p_vec[:,0]==0]
        # ww_canti = ww_canti[len(ww_canti)//2:len(ww_canti)]
        # I_canti = I_vec[p_vec[:,0]==0]
        # I_canti = I_canti[len(I_canti)//2:len(I_canti)]
        
        y_canti = y[len(y)//2:len(y)]
        b_canti = max(y_canti) - min(y_canti)
        l_canti = l_vec[len(y)//2:len(y)]*dy
        ww_canti = (p_vec[:,1]+ww_vec[:,0])[len(y)//2:len(y)]
        I_canti = I_vec[len(y)//2:len(y)]
        
        L_canti = sum(l_canti)
        # Generate Local Stiffness Matrices
        n_canti = len(l_canti)
        k_l = np.zeros((4,4,n_canti))
        for i in range(n_canti):
            k_l[:,:,i] = (E*I_canti[i]/dy**3)*np.array([[  12,  6*dy   ,   -12, 6*dy   ],
                                                        [6*dy,  4*dy**2, -6*dy, 2*dy**2],
                                                        [ -12, -6*dy   , 12   ,-6*dy   ],
                                                        [6*dy,  2*dy**2, -6*dy, 4*dy**2]])
        
        # Assemble Global Stiffness Matrix
        k_g = np.zeros((2*n_canti,2*n_canti))
        for i in range(n_canti):
            if i == 0:
                k_g[0:2,0:2] = k_l[2:4,2:4,i]
            else:
                k_g[2*(i-1):2*(i+1),2*(i-1):2*(i+1)] = k_g[2*(i-1):2*(i+1),2*(i-1):2*(i+1)] + k_l[:,:,i]
        
        k_g_inv = np.linalg.inv(k_g)
        
        # Assemble Force Vectors
        tipLoadOffset = 1.12*25.4/1000
        force_vec = np.zeros((2*n_canti,3))
        midPointNotFound = 1
        tipLoadPointNotFound = 1
        for i in range(n_canti):
            force_vec[2*i,0] = l_canti[i]
            if midPointNotFound and y_canti[i] <= b/4 and y_canti[i+1] > b/4:
                force_vec[2*i,1] = wLoad
                midPointNotFound = 0
                print("MIDPOINT FOUND TO APPLY LOAD!!!!!!!!")
            if tipLoadPointNotFound and y_canti[i] >= b/2 - tipLoadOffset:
                force_vec[2*i,2] = wLoad
                tipLoadPointNotFound = 0
                print("TIP LOAD POINT FOUND!!!!!!!!")
            # if i == n_canti-1:
            #     force_vec[2*i,2] = wLoad
            
        disp_vec = np.zeros((2*n_canti,8))
        
        disp_vec[:,0] = k_g_inv@force_vec[:,0]      # Lift Dist 1g
        disp_vec[:,1] = k_g_inv@(3*force_vec[:,0])  # Lift Dist 3g
        disp_vec[:,2] = k_g_inv@(5*force_vec[:,0])  # Lift Dist 5g
    
        disp_vec[:,3] = k_g_inv@force_vec[:,1]                   # Mid-Load 1 Bag
        disp_vec[:,4] = k_g_inv@(2*force_vec[:,1])               # Mid-Load 2 Bags
        disp_vec[:,5] = k_g_inv@force_vec[:,2]                   # Tip-Load 1 Bag
        disp_vec[:,6] = k_g_inv@(force_vec[:,1] + force_vec[:,2])   # 1 Mid, 1 Tip 
        disp_vec[:,7] = k_g_inv@(2*force_vec[:,2])               # Tip-Load 2 Bags
        
        
        
        delta_lin = np.zeros((n_canti,8))
        delta_ang = np.zeros((n_canti,8))
        for i in range(n_canti):
            # print()
            delta_lin[i,:] = disp_vec[2*i,:]
            delta_ang[i,:] = disp_vec[2*i+1,:]
        
        
        deltaTip = 0.5*L*(b_canti)**3/(3*E*I_vec[n//2]) # Tip Deflection, loaded at tip
        deltaMid = (((0.5*L)*(b_canti/2)**2)/(6*E*I_vec[n//2]))*(3*b_canti - b_canti/2) # Tip Deflection, loaded at mid point
        
        m2in = 39.37
        m2mm = 1000
        
        print("**Bend Predictions**")
        print("Mid Load 1 Bag Analytic Def = %5.1f mm = %4.2f in" %(deltaMid*m2mm,deltaMid*m2in))
        print("Mid Load 1 Bag FEM Def      = %5.1f mm = %4.2f in" %(delta_lin[-1,3]*m2mm,delta_lin[-1,3]*m2in))
        print("Mid Load 2 Bags FEM Def     = %5.1f mm = %4.2f in" %(delta_lin[-1,4]*m2mm,delta_lin[-1,4]*m2in))
        print("Tip Load 1 Bag Analytic Def = %5.1f mm = %4.2f in" %(deltaTip*m2mm,deltaTip*m2in))
        print("Tip Load 1 Bag FEM Def      = %5.1f mm = %4.2f in" %(delta_lin[-1,5]*m2mm,delta_lin[-1,5]*m2in))
        print("Mid Bag + Tip Bag FEM Def   = %5.1f mm = %4.2f in" %(delta_lin[-1,6]*m2mm,delta_lin[-1,6]*m2in))
        print("Tip Load 2 Bags FEM Def     = %5.1f mm = %4.2f in" %(delta_lin[-1,7]*m2mm,delta_lin[-1,7]*m2in))
        print("********************")
        
        bendPredict = open("WingBend.csv", "a")  # append mode
        for i in range(len(y_canti)):
            if i!=0:
                bendPredict.write("\n")
            bendPredict.write("%5.4e," %y_canti[i])
            bendPredict.write("%5.4e," %delta_lin[i,3])
            bendPredict.write("%5.4e," %delta_lin[i,4])
            bendPredict.write("%5.4e," %delta_lin[i,5])
            bendPredict.write("%5.4e," %delta_lin[i,6])
            bendPredict.write("%5.4e" %delta_lin[i,7])
        bendPredict.close()
        
        if dispPlots:
        
            width  = 6.5
            height = 5
            
            # plt.subplot(212)
            plt.figure(figsize=(width, height))
            plt.plot(y_canti*1000,1000*delta_lin[:,0],'-',color = 'k',label = r'FEM Lift Dist @ 1g')
            plt.plot(y_canti*1000,1000*delta_lin[:,1],':',color = 'k',label = r'FEM Lift Dist @ 3g')
            plt.plot(y_canti*1000,1000*delta_lin[:,2],'--',color = 'k',label = r'FEM Lift Dist @ 5g')
            plt.plot(0.5*b*1000,deltaMid*1000,'*',color = 'k',label = r'An. Mid Load @ 1g')
            plt.grid()
            plt.ylabel(r'Vertical Deflection, mm', fontsize=20)
            # plt.xlim(0,0.5*b*1000)
            
            plt.xticks(fontsize = 18)
            plt.yticks(fontsize = 18)
            plt.xlabel(r'Spanwise Location, mm', fontsize=20)
            plt.ylim(None,1.25*max(1000*delta_lin[:,2]))
            # plt.legend(fontsize=9)
            plt.legend(fontsize=16,bbox_to_anchor=(0., 0.9, 0.65, .10), loc=2,
                ncol=1, mode="expand", borderaxespad=0.)
            
            width  = 6.5
            height = 5
            
            # plt.subplot(212)
            plt.figure(figsize=(width, height))
            plt.plot(y_canti*1000,1000*delta_lin[:,3],'--',color = 'black',   label = r'FEM Mid Load 1 Bag')
            plt.plot(y_canti*1000,1000*delta_lin[:,4],'--',color = 'blue',label = r'FEM Mid Load 2 Bags')
            plt.plot(y_canti*1000,1000*delta_lin[:,5],'--',color = 'green',   label = r'FEM Tip Load 1 Bag')
            plt.plot(y_canti*1000,1000*delta_lin[:,6],'--',color = 'orange',label = r'FEM 1 Mid + 1 Tip')
            plt.plot(y_canti*1000,1000*delta_lin[:,7],'--',color = 'red',label = r'FEM Tip Load 2 Bags')
            plt.plot(0.5*b*1000,deltaMid*1000,'*',color = 'black',label = r'An. Mid Load 1 Bag')
            plt.plot(0.5*b*1000,deltaTip*1000,'*',color = 'green',label = r'An. Tip Load 1 Bag')
            plt.grid()
            plt.ylabel(r'Vertical Deflection, mm', fontsize=20)
            # plt.xlim(0,0.5*b*1000)
            
            plt.xticks(fontsize = 18)
            plt.yticks(fontsize = 18)
            plt.xlabel(r'Spanwise Location, mm', fontsize=20)
            plt.ylim(None,1.5*max(1000*delta_lin[:,2]))
            # plt.legend(fontsize=9)
            plt.legend(fontsize=16,bbox_to_anchor=(0., 0.9, 0.65, .10), loc=2,
                ncol=1, mode="expand", borderaxespad=0.)
            
            
            plt.figure(figsize=(width, height))
            plt.plot(y_canti*39.37,1000*delta_lin[:,3]/25.4,'--',color = 'black',   label = r'FEM Mid Load 1 Bag')
            plt.plot(y_canti*1000/25.4,1000*delta_lin[:,4]/25.4,'--',color = 'blue',label = r'FEM Mid Load 2 Bags')
            plt.plot(y_canti*1000/25.4,1000*delta_lin[:,5]/25.4,'--',color = 'green',   label = r'FEM Tip Load 1 Bag')
            plt.plot(y_canti*1000/25.4,1000*delta_lin[:,6]/25.4,'--',color = 'orange',label = r'FEM 1 Mid + 1 Tip')
            plt.plot(y_canti*1000/25.4,1000*delta_lin[:,7]/25.4,'--',color = 'red',label = r'FEM Tip Load 2 Bags')
            plt.plot(0.5*b*1000/25.4,deltaMid*1000/25.4,'*',color = 'black',label = r'An. Mid Load 1 Bag')
            plt.plot(0.5*b*1000/25.4,deltaTip*1000/25.4,'*',color = 'green',label = r'An. Tip Load 1 Bag')
            plt.grid()
            plt.ylabel(r'Vertical Deflection, in', fontsize=20)
            # plt.xlim(0,0.5*b*1000)
            
            plt.xticks(fontsize = 18)
            plt.yticks(fontsize = 18)
            plt.xlabel(r'Spanwise Location, in', fontsize=20)
            plt.ylim(None,max(1000*delta_lin[:,7]/25.4))
            plt.legend(fontsize=16,bbox_to_anchor=(0., 0.9, 0.65, .10), loc=2,
                ncol=1, mode="expand", borderaxespad=0.)
        fout = [delta_lin[-1,0],delta_lin[-1,2],delta_lin[-1,3]]
    except:
        print("Error in Wing Bend Calculation... Skipped!")
        fout = [np.nan,np.nan,np.nan]
        
    return(fout)