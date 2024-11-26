def ChordTwist(n,b,AR,cL1,cL2,a,rho,W_G,planform,twist,dispPlots):
    """
    This function allows the user to directly specify chord and twist
        distributions.  The function will then use LL-theory to calculate the
        AoA needed for cruise (cL1) and pitch and wing to that AoA.  The
        function also used LL-theory to calculate the lift distribution at 
    
    Parameters
    ----------
    n : Scalar
        Number of discretized spanwise points
    b : Scalar
        Wing Span
    AR : Scalar
        Aspect Ratio of Wing
    cL1 : Scalar
        cL in cruise, the wing alpha dist will be returned for this cL
    cL2 : Scalar
        cL max, for lift distribution for structural design
    a : Scalar
        airfoil lift curve slope
    rho : Scalar
        Air density
    W_G : Scalar
        Gross Aircraft Weight
    planform : 3 possible types
        'e'     can be given for an elliptical planform
        scalar  can be given to define the taper ratio of a single taper wing
        mx2 numpy array can be given to specify taper breaks at m wing points
            if m[0,0] = 0: wing is symmetric, specifying right wing only
            if m[0,0] =-0.5: entire span is defined 
    twist : mx2 numpy array
        scalar  can be given to define the washout for a linear twist
        mx2 numpy array can be given to specify twist breaks at m wing points
            if m[0,0] = 0: wing is symmetric, specifying right wing only
            if m[0,0] =-0.5: entire span is defined 

    Returns
    -------
    WingShape : nx3 numpy array
        Columns: (0) spanwise location (y/b)
                 (1) chord (c/b)
                 (2) alpha (rad) (@cL1, used for mounting angle)
    l_vec : nx1 numpy array
        lift distribution at cL2 used for structural design

    """
    import numpy as np
    from numpy.linalg import inv
    from numpy import pi
    
    # b = b/1000 #convert mm in m
    
    # Define needed vectors and matrices ------------------------
    e = np.ones((n,1))
    I = np.eye(n)
    
    yF = np.linspace(0.5*(-1 + 1/n), 0.5*(1 - 1/n), n)
    y = np.zeros((n,1))
    for i in range(n):
        y[i,0] = yF[i]
    del(yF)
    Q = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            Q[i,j] = 1/(1-4*(i-j)**2)
    Qinv = inv(Q)
    print("Q-Matrix Generated!\n")
    W = ((a*n)/(2*pi))*Q
    
    # Define Chord Shape ----------------------------------
    if planform == 'e': # Elliptical Wing
        c_sh = Qinv@e
    elif isinstance(planform, float) or isinstance(planform, int): # Single Taper
        singleWingShape = np.linspace(planform, 1, n//2)
        c_sh = np.concatenate((singleWingShape,np.flip(singleWingShape)))
    else:
        print("CHORD INPUT ERROR!!!!")
        import sys
        sys.exit()
    
    # Convert to column Vector
    c_shTemp = np.zeros((n,1))
    for i in range(n):
        c_shTemp[i,0] = c_sh[i]
    c_sh = c_shTemp
    del(c_shTemp)

    c_vec = (n/AR)*(c_sh/(e.T@c_sh))
    C     = np.diag(np.ndarray.flatten(c_vec))
    Cinv  = np.diag(np.ndarray.flatten(1/c_vec))
    
    if isinstance(twist, float) or isinstance(twist, int): # Linear Washout/in
        rightWingTwist = np.linspace(0, twist, n//2)
        leftWingTwist  = np.flip(rightWingTwist)
        alphatw_vec = (pi/180)*np.concatenate((leftWingTwist,rightWingTwist))
    else:
        print("TWIST INPUT ERROR!")
        import sys
        sys.exit()
        
        
    # Convert to column Vector
    alphatwTemp = np.zeros((n,1))
    for i in range(n):
        alphatwTemp[i,0] = alphatw_vec[i]
    alphatw_vec = alphatwTemp
    del(alphatwTemp)
    
    # Lifting Line Math
    cLalpha   = ((AR/n)*e.T@inv(Cinv + W)@e*a)[0,0] #3D lift curve slope
    cL0       = (a*AR/n)*(e.T@inv(Cinv+W)@alphatw_vec)[0,0]   #cL @ alpha_w = 0
    
    # Calculations for cruise
    V1           = np.sqrt(2*W_G/(rho*(b**2/AR)*cL1))   
    alpha_w1     = (1/cLalpha)*(cL1 - cL0)
    Gamma_vec1   = 0.5*a*inv(Cinv + W)@(alphatw_vec + alpha_w1)
    l_vec1       = rho*V1*Gamma_vec1
    alpha_i_vec1 = (n/pi)*Q@Gamma_vec1
    
    # Calculations for cLmax    
    V2         = np.sqrt(2*W_G/(rho*(b**2/AR)*cL2))
    alpha_w2   = (1/cLalpha)*(cL2 - cL0)
    Gamma_vec2 = 0.5*V2*b*a*inv(Cinv + W)@(alphatw_vec + alpha_w2)
    l_vec2     = rho*V2*Gamma_vec2
        
    OutputString ="""|------------------------------------------|
|  cLalpha = %3.2f*a   ;    cL0 = %3.2f     |
|------------------------------------------|
| Condition  |   Cruise    |    cLmax      |
| cL         |   %4.3f     |    %4.3f      |
| alpha_w    |   %4.3f deg |    %4.3f deg |
|------------------------------------------|\n""" \
    % (cLalpha/a,cL0,cL1, cL2, (180/pi)*alpha_w1 , (180/pi)*alpha_w2 )
    print(OutputString)
    del(OutputString)

    # Plotting --------------------------------------------------
    if dispPlots == 1:
        import matplotlib.pyplot as plt
        
        width  = 6.5
        height = 5.5
        plt.figure(figsize=(width, height))
        
        plt.subplot(211)
        plt.grid()
        plt.ylabel(r'Chord, $c/b$', fontsize=20)
        # plt.axis('equal')
        plt.plot(y, c_vec,'k-')
        plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)), (), fontsize = 18)
        plt.yticks(np.array((0,0.05,0.10)), fontsize = 18)
        plt.ylim(0, 1.2*max(c_vec))
        plt.xlim(-0.5,0.5)
        
        plt.subplot(212)
        plt.grid()
        plt.ylabel('Angle of\n'+' Attack, deg', fontsize=20)
        plt.ylim(-0.2, 1.5*(180/pi)*max((alphatw_vec + alpha_w1)))
        plt.xlim(-0.5,0.5)
        plt.plot(y, (180/pi)*alphatw_vec,'k--', label = r'$\alpha_{tw}$')
        plt.plot(y, (180/pi)*(alphatw_vec + alpha_w1),'k-', label = r'$\alpha_{tw} + \alpha_w$')
        plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)), fontsize = 18)
        # plt.yticks(np.array((0,-1,-2)), fontsize = 18)
        plt.yticks(fontsize = 18)
        plt.legend(fontsize=18,bbox_to_anchor=(0., 0.9, 1., .10), loc=2,
           ncol=3, mode="expand", borderaxespad=0.)
        plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
        
        width  = 6.5
        height = 5.5
        plt.figure(figsize=(width, height))
        
        plt.subplot(211)
        plt.grid()
        plt.ylabel(r'Vorticity, $\frac{\Gamma}{V_\infty b}$', fontsize=20)
        # plt.ylim(bottom=0)
        plt.xlim(-0.5,0.5)
        plt.plot(y, Gamma_vec1,'k-')
        # plt.plot(y, Gamma_vec2,'k-')
        plt.yticks(fontsize = 18)
        plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)), (), fontsize = 18)
        
        plt.subplot(212)
        plt.grid()
        plt.ylabel(r'Downwash, $\frac{w}{V_\infty}$', fontsize=20)
        # plt.ylim(bottom=0)
        plt.xlim(-0.5,0.5)
        plt.plot(y, alpha_i_vec1,'k-')
        # plt.plot(y, Gamma_vec2,'k-')
        plt.yticks(fontsize = 18)
        plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)), fontsize = 18)
        plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
        
        # plt.subplot(414)
        # plt.grid()
        # plt.ylabel('Normalized\nLift')
        # # plt.ylim(bottom=0)
        # plt.plot(y, l_vec1*b/W_G,'k--',label = r'c_L = '+str(cL1))
        # plt.plot(y, l_vec2*b/W_G,'k-', label = r'c_L = '+str(cL2))
        # plt.xlim(-0.5,0.5)
        # plt.xlabel(r'Spanwise Location, $y/b$', fontsize=10)
        # plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)))
        # plt.legend(fontsize=9)
    
      
    # Build Wing Shape Array
    WingShape = np.concatenate((y,c_vec,alphatw_vec+alpha_w1),axis=1) 
    
    return [WingShape,l_vec2 ]


def LiftChord(n,b,AR,cL1,cL2,a,rho,W_G,liftShape,planform,dispPlots):
    """
    This function allows the user to specify lift and chord distributions.
    
    Parameters
    ----------
    n : Scalar
        Number of discretized spanwise points
    b : Scalar
        Wing Span
    AR : Scalar
        Aspect Ratio of Wing
    cL1 : Scalar
        cL in cruise, the wing alpha dist will be returned for this cL
    cL2 : Scalar
        cL max, for lift distribution for structural design
    a : Scalar
        airfoil lift curve slope
    rho : Scalar
        Air density
    W_G : Scalar
        Gross Aircraft Weight
    planform : 3 possible types
        'e'     can be given for an elliptical planform
        scalar  can be given to define the taper ratio of a single taper wing
        mx2 numpy array can be given to specify taper breaks at m wing points
            if m[0,0] = 0: wing is symmetric, specifying right wing only
            if m[0,0] =-0.5: entire span is defined 
    twist : mx2 numpy array
        scalar  can be given to define the washout for a linear twist
        mx2 numpy array can be given to specify twist breaks at m wing points
            if m[0,0] = 0: wing is symmetric, specifying right wing only
            if m[0,0] =-0.5: entire span is defined 

    Returns
    -------
    WingShape : nx3 numpy array
        Columns: (0) spanwise location (y/b)
                 (1) chord (c/b)
                 (2) alpha (rad) (@cL1, used for mounting angle)
    l_vec2 : nx1 numpy array
        lift distribution at cL2 used for structural design

    """
    import numpy as np
    from numpy.linalg import inv
    from numpy import pi
    # b = b/1000
    
    # Define needed vectors and matrices ------------------------
    e = np.ones((n,1))
    I = np.eye(n)
    
    yF = np.linspace(0.5*(-1 + 1/n), 0.5*(1 - 1/n), n)
    y = np.zeros((n,1))
    for i in range(n):
        y[i,0] = yF[i]
    del(yF)
    
    yF = np.linspace(0.5*(-1 + 1/n), 0.5*(1 - 1/n), n)
    y = np.zeros((n,1))
    for i in range(n):
        y[i,0] = yF[i]
    del(yF)
    Q = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            Q[i,j] = 1/(1-4*(i-j)**2)
    Qinv = inv(Q)
    print("Q-Matrix Generated!\n")
    W = ((a*n)/(2*pi))*Q
    
    A = a*I
    Ainv = (1/a)*I
    
    # Define Aero Design Lift Distribution ----------------------------------
    if liftShape == 'e':
        l_sh = Qinv@e
        
    V1         = np.sqrt(2*W_G/(rho*(b**2/AR)*cL1))   
    Gamma_vec1 = (1/2)*(cL1/AR)*(n*l_sh/(e.T@l_sh))
    l_vec1     = rho*V1*Gamma_vec1
    
    # Define Chord Shape ----------------------------------
    if planform == 'e': # Elliptical Wing
        c_sh = Qinv@e
    elif isinstance(planform, float) or isinstance(planform, int): # Single Taper
        singleWingShape = np.linspace(planform, 1, n//2)
        c_sh = np.concatenate((singleWingShape,np.flip(singleWingShape)))
    elif isinstance(planform, list) and planform[0] == 'skewedEllipse':
        skewArray = np.array([np.linspace(1-0.5*planform[1],1+0.5*planform[1],n)]).T
        c_sh = (Qinv@e)*skewArray
    elif isinstance(planform, list) and planform[0] == 'cuffed':
        # planform = ['cuffed', TaperRatio, yCuff, ChordDelta]
        TR  = planform[1]
        y_c = planform[2]
        d_c = planform[3]
        c_r = 1
        y_t = 0.5
        c_t = c_r*TR
        c_c = ((y_c/y_t)*(c_t - c_r) + c_r)/(1+(y_c/y_t)*d_c)
        chordGrad = ((c_t - d_c - c_r)/y_t)
        
        c_sh = np.zeros((n,1))
        for j in range(n):
            if abs(y[j,0]) > y_c*y_t:
                isTip = 1
            else:
                isTip = 0
            c_sh[j,0] = np.sign(y[j,0])*chordGrad*y[j,0] + c_r + isTip*d_c*c_c
    else:
        print("CHORD INPUT ERROR!!!!")
        import sys
        sys.exit()
    
    # Convert to column Vector
    c_shTemp = np.zeros((n,1))
    for i in range(n):
        c_shTemp[i,0] = c_sh[i]
    c_sh = c_shTemp
    del(c_shTemp)
    c_vec      = (n/AR)*(c_sh/(e.T@c_sh))
    Cinv = np.diag(1/np.ndarray.flatten(c_vec))    
    
    alpha_vec1 = 2*Ainv@(Cinv + W)@Gamma_vec1
    alpha_w1 = alpha_vec1[n//2,0]
    alphatw_vec = alpha_vec1 - alpha_w1
    
    # Lifting Line Math
    cLalpha    =  a*((e.T@inv(Cinv + W)@e)/(e.T@c_vec))[0,0] #3D lift curve slope
    cL0        = (a*AR/n)*(e.T@inv(Cinv+W)@alphatw_vec)[0,0]   #cL @ alpha_w = 0
    alpha_w0   = alpha_w1 - (cL1/cLalpha)
    alpha_vec0 = alphatw_vec + alpha_w0
    Gamma_vec0 = (1/2)*inv(Cinv + W)@A@alpha_vec0
    
    V2         = np.sqrt(2*W_G/(rho*(b**2/AR)*cL2))
    alpha_w2   = (1/cLalpha)*(cL2-cL1) + alpha_w1
    alpha_vec2 = alphatw_vec + alpha_w2
    Gamma_vec2 = (1/2)*inv(Cinv + W)@A@alpha_vec2
    l_vec2     = rho*V2*Gamma_vec2
    
    
    
    # Gamma_vec2 = 0.5*V2*b*a*inv(Cinv + W)@(alphatw_vec + alpha_w2)
    # breakpoint()
    
    OutputString ="""|------------------------------------------|
|  cLalpha = %3.2f*a   ;    cL0 = %3.2f     |
|------------------------------------------|
| Condition  |   Cruise    |    cLmax      |
| cL         |   %4.3f     |    %4.3f      |
| V          |   %4.2f m/s |    %4.2f m/s  |
| alpha_w    |   %4.3f deg |    %4.3f deg |
|------------------------------------------|\n""" \
    % (cLalpha/a,cL0,cL1, cL2, V1,V2,(180/pi)*alpha_w1 , (180/pi)*alpha_w2 )
    print(OutputString)
    del(OutputString)
    if dispPlots:
        # Plotting --------------------------------------------------
        import matplotlib.pyplot as plt
        
        width  = 6.5
        height = 5.5
        plt.figure(figsize=(width, height))
        
        plt.subplot(211)
        plt.grid()
        plt.ylabel(r'Chord, $c/b$', fontsize=20)
        # plt.axis('equal')
        plt.plot(y, c_vec,'k-')
        plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)), (), fontsize=18)
        plt.yticks(np.array((0,0.05,0.10)), fontsize=18)
        plt.ylim(0, 1.2*max(c_vec))
        plt.xlim(-0.5,0.5)
        
        plt.subplot(212)
        plt.grid()
        plt.ylabel(r'Vorticity, $\frac{\Gamma}{V_\infty b}$', fontsize=20)
        # plt.ylim(bottom=0)
        plt.xlim(-0.5,0.5)
        plt.plot(y, Gamma_vec1,'k-')
        # plt.plot(y, Gamma_vec2,'k-')
        plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)), fontsize=18)
        plt.yticks(fontsize=18)
        plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
        
        
        width  = 6.5
        height = 5.5
        plt.figure(figsize=(width, height))
        
        plt.subplot(212)
        plt.grid()
        plt.ylabel('Angle of\n'+' Attack, deg', fontsize=20)
        # plt.ylim(0, 1.2*(180/pi)*max(alpha))
        plt.xlim(-0.5,0.5)
        plt.plot(y, (180/pi)*alphatw_vec,'k--', label=r'$\alpha_\mathrm{tw}$')
        plt.plot(y, (180/pi)*alpha_vec1,'k-', label=r'$\alpha_\mathrm{tw}+\alpha_w$')
        plt.yticks(fontsize=18)
        plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
        plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)), fontsize=18)
        plt.legend(fontsize=18,ncol = 2)
        
        
        width  = 6.5
        height = 5.5
        plt.figure(figsize=(width, height))
        plt.subplot(211)
        plt.grid()
        plt.ylabel('Angle of\n'+' Attack, deg', fontsize=20)
        # plt.ylim((180/pi)*min(min(alphatw_vec),min(alpha_vec1),min(alpha_vec2)), 1.5*(180/pi)*max(max(alphatw_vec),max(alpha_vec1),max(alpha_vec2)))
        plt.xlim(-0.5,0.5)
        plt.plot(y, (180/pi)*alpha_vec2,'k:', label=r'$c_L =$ %3.2f'%(cL2))
        plt.plot(y, (180/pi)*alpha_vec1,'k-', label=r'$c_L =$ %3.2f'%(cL1))
        plt.plot(y, (180/pi)*alpha_vec0,'k--', label=r'$c_L = 0$')
        plt.yticks(fontsize=18)
        # plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
        plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)), (), fontsize=16)
        plt.legend(fontsize=18,ncol = 1,bbox_to_anchor=(1.0, .3))
        
        plt.subplot(212)
        plt.grid()
        plt.ylabel(r'Vorticity, $\frac{\Gamma}{V_\infty b}$', fontsize=20)
        # plt.ylim(bottom=0)
        plt.xlim(-0.5,0.5)
        plt.plot(y, Gamma_vec0,'k--')
        plt.plot(y, Gamma_vec1,'k-')
        plt.plot(y, Gamma_vec2,'k:')
        # plt.plot(y, Gamma_vec2,'k-')
        plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)), fontsize=18)
        plt.yticks(fontsize=18)
        plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
        
        
        
        # plt.subplot(414)
        # plt.grid()
        # plt.ylabel('Normalized\nLift')
        # # plt.ylim(bottom=0)
        # plt.plot(y, l_vec1*b/W_G,'k--',label = r'c_L = '+str(cL1))
        # plt.plot(y, l_vec2*b/W_G,'k-', label = r'c_L = '+str(cL2))
        # plt.xlim(-0.5,0.5)
        # plt.xlabel(r'Spanwise Location, $y/b$', fontsize=10)
        # plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)))
        # plt.legend(fontsize=9)
    
    # breakpoint()
    
    # Return Wing Shape
    WingShape = np.concatenate((y,c_vec,alpha_vec1),axis=1)    
    return [WingShape, l_vec2]


def LiftTwist(n,b,AR,cL1,cL2,a,rho,W_G,liftShape,twist,fileName,nAVL,xcBV,dispPlots):
    """
    This function allows the user to specify lift and chord distributions.
    
    Parameters
    ----------
    n : Scalar
        Number of discretized spanwise points
    b : Scalar
        Wing Span (m)
    AR : Scalar
        Aspect Ratio of Wing
    cL1 : Scalar
        cL in cruise, the wing alpha dist will be returned for this cL
    cL2 : Scalar
        cL max, for lift distribution for structural design
    a : Scalar
        airfoil lift curve slope (cl/rad)
    rho : Scalar
        Air density 
    W_G : Scalar
        Gross Aircraft Weight 
    planform : 3 possible types
        'e'     can be given for an elliptical planform
        scalar  can be given to define the taper ratio of a single taper wing
        mx2 numpy array can be given to specify taper breaks at m wing points
            if m[0,0] = 0: wing is symmetric, specifying right wing only
            if m[0,0] =-0.5: entire span is defined 
    twist : mx2 numpy array
        scalar  can be given to define the washout for a linear twist
        mx2 numpy array can be given to specify twist breaks at m wing points
            if m[0,0] = 0: wing is symmetric, specifying right wing only
            if m[0,0] =-0.5: entire span is defined 

    Returns
    -------
    WingShape : nx3 numpy array
        Columns: (0) spanwise location (y/b)
                 (1) chord (c/b)
                 (2) alpha (rad) (@cL1, used for mounting angle)
    l_vec2 : nx1 numpy array
        lift distribution at cL2 used for structural design

    """
    import numpy as np
    from numpy.linalg import inv
    from numpy import pi
    import copy
    
    # Define needed vectors and matrices ------------------------
    e = np.ones((n,1))
    I = np.eye(n)
    
    yF = np.linspace(0.5*(-1 + 1/n), 0.5*(1 - 1/n), n)
    y = np.zeros((n,1))
    for i in range(n):
        y[i,0] = yF[i]
    del(yF)
    
    yF = np.linspace(0.5*(-1 + 1/n), 0.5*(1 - 1/n), n)
    y = np.zeros((n,1))
    for i in range(n):
        y[i,0] = yF[i]
    del(yF)
    Q = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            Q[i,j] = 1/(1-4*(i-j)**2)
    Qinv = inv(Q)
    print("Q-Matrix Generated!\n")
    W = ((a*n)/(2*pi))*Q
    
    A = a*I
    Ainv = (1/a)*I
    
    # Define Aero Design Lift Distribution ----------------------------------
    if liftShape == 'e':
        l_sh = Qinv@e
    elif isinstance(liftShape, list) and liftShape[0] == 'b':
        eta = liftShape[1]
        K = np.zeros((n,n))
        for j in range(n):
            for k in range(n):
                K[j,k] = (1/n**2)*abs(j-k)
        l_sh = Qinv@(I - eta*K)@e
        
    V1         = np.sqrt(2*W_G/(rho*(b**2/AR)*cL1))   
    Gamma_vec1 = (1/2)*(n*cL1/AR)*(l_sh/(e.T@l_sh))
    l_vec1     = rho*V1*Gamma_vec1
    
    # Define Twist Shape ----------------------------------
    if twist == 'e': # Elliptical Wing
        tw_sh = Qinv@e
    elif isinstance(twist, float) or isinstance(twist, int): # Single Taper
        singleTwistShape = np.linspace(twist*np.pi/180, 0, n//2)
        tw_sh = np.concatenate((singleTwistShape,np.flip(singleTwistShape)))
    elif isinstance(twist, list) and twist[0] == 'assym':
        tw_sh = (pi/180)*np.linspace(-twist[1], twist[1], n)
    elif isinstance(twist, list) and twist[0] == 'skewedEllipse':
        skewArray = (pi/180)*np.array([np.linspace(1-0.5*twist[1],1+0.5*twist[1],n)]).T
        ell_vec = Qinv@e
        ell_vec_norm = n*ell_vec/(e.T@ell_vec)
        tw_sh = ell_vec_norm*skewArray
    elif isinstance(twist, list) and twist[0] == 'sine':
        tw_sh = (pi/180)*twist[1]*np.sin(2*pi*y)
    else:
        print("TWIST INPUT ERROR!!!!")
        import sys
        sys.exit()
    
    # Convert to column Vector
    tw_shTemp = np.zeros((n,1))
    for i in range(n):
        tw_shTemp[i,0] = tw_sh[i]
    alphatw_vec = tw_shTemp
    del(tw_shTemp)
    
    limTerm = 2*Ainv@W@Gamma_vec1 - alphatw_vec
    if any(Gamma_vec1>0) and any(Gamma_vec1<0):
        print("Lift Pos and Neg")
        Jpos = Gamma_vec1 >= 0
        Jneg = Gamma_vec1 < 0
        alpha_w_min = max(limTerm[Jpos])
        alpha_w_max = min(limTerm[Jneg])
        print("alpha_w_min = %5.3f deg" %((180/pi)*alpha_w_min))
        print("alpha_w_max = %5.3f deg" %((180/pi)*alpha_w_max))
        
        # breakpoint()
        
        
    elif all(Gamma_vec1>0):
        print("Lift All Pos")
        alpha_w_min = max(limTerm)[0]
        alpha_w_max = alpha_w_min + 2*AR*max(Ainv@Gamma_vec1)[0]
        
    elif all(Gamma_vec1<0):
        print("Lift All Neg")
        alpha_w_max = min(limTerm)
        alpha_w_min = alpha_w_max - 2*AR*min(Ainv@Gamma_vec1)
        
    else:
        print("ALPHA_W LIMIT ERROR!!!!")
        import sys
        sys.exit()
    alpha_w_Lims0 = np.array([alpha_w_min, alpha_w_max])
    dalphaLims = alpha_w_max - alpha_w_min
    
    alpha_w_sweepLims = np.array([alpha_w_min - dalphaLims,alpha_w_max + dalphaLims])
    alpha_w_sweep = np.linspace(alpha_w_sweepLims[0],alpha_w_sweepLims[1],1001)
    AR_sweep = np.zeros(len(alpha_w_sweep))
    count = -1
    for alpha_w in alpha_w_sweep:
        count += 1
        alpha_vec = alphatw_vec + alpha_w
        diagMat = np.diag(1/np.ndarray.flatten(0.5*A@alpha_vec - W@Gamma_vec1))
        AR_sweep[count] = n/(e.T@diagMat@Gamma_vec1)
    
    
    # alpha_w_Lims0 = np.array([0,10])*(np.pi/180)
    alpha_w_Lims = copy.deepcopy(alpha_w_Lims0)
    alphaTol = 1e-4
    alphaNotFound=1
    count = 0
    while alphaNotFound and count < 100:
        count += 1
        print("i = ", count)
        
        alpha_w = 0.5*(alpha_w_Lims[0] + alpha_w_Lims[1])
        alpha_vec = alphatw_vec + alpha_w
    #     print(alpha_vec)
    
        # Calculate chord
        Aalpha = A@alpha_vec
        Wl     = W@Gamma_vec1
        c_vec = np.zeros((n,1))
        for j in range(n):
            c_vec[j,0] = Gamma_vec1[j,0]/(0.5*Aalpha[j,0] - Wl[j,0])
    #     print(c_vec)
        
        ARcalc = (n/(e.T@c_vec))[0,0]
        ARres  = (AR - ARcalc)/AR
    #     print("alpha_w = ", alpha_w*180/np.pi,"AR_calc = ",ARcalc, ";     AR_res = ",ARres)
        
        if abs(ARres) < alphaTol:
            alphaNotFound = 0
        else:
            if ARres > 0:
                alpha_w_Lims[0] = alpha_w
            else:
                alpha_w_Lims[1] = alpha_w
    
    print("Count = ", count, ";   alpha_w = ", alpha_w*180/np.pi)
    print("AR_calc = ",ARcalc, ";     AR_res = ",ARres)

    alpha_w1   = alpha_w
    alpha_vec1 = alpha_vec
    
    
    # c_vec      = (n/AR)*(c_sh/(e.T@c_sh))
    Cinv = np.diag(1/np.ndarray.flatten(c_vec))    
    
    # alpha_vec1 = 2*Ainv@(Cinv + W)@Gamma_vec1
    # alpha_w1 = alpha_vec1[n//2,0]
    # alphatw_vec = alpha_vec1 - alpha_w1
    
    # Lifting Line Math
    cLalpha    =  a*((e.T@inv(Cinv + W)@e)/(e.T@c_vec))[0,0] #3D lift curve slope
    cL0        = (a*AR/n)*(e.T@inv(Cinv+W)@alphatw_vec)[0,0]   #cL @ alpha_w = 0
    V2         = np.sqrt(2*W_G/(rho*(b**2/AR)*cL2))
    alpha_w2   = (1/cLalpha)*(cL2-cL1) + alpha_w1
    alpha_vec2 = alphatw_vec + alpha_w2
    Gamma_vec2 = (1/2)*inv(Cinv + W)@A@alpha_vec2
    l_vec2     = rho*V2*Gamma_vec2
    
    Gamma_vec2 = 0.5*V2*b*a*inv(Cinv + W)@(alphatw_vec + alpha_w2)
    # breakpoint()
    
    # Return Wing Shape
    WingShape = np.concatenate((y,c_vec,alpha_vec1),axis=1)    

    
    
    
    OutputString ="""|------------------------------------------|
|  cLalpha = %3.2f*a   ;    cL0 = %3.2f     |
|------------------------------------------|
| Condition  |   Cruise    |    cLmax      |
| cL         |   %4.3f     |    %4.3f      |
| V          |   %4.2f m/s |    %4.2f m/s  |
| alpha_w    |   %4.3f deg |    %4.3f deg |
|------------------------------------------|\n""" \
    % (cLalpha/a,cL0,cL1, cL2, V1,V2,(180/pi)*alpha_w1 , (180/pi)*alpha_w2 )
    print(OutputString)
    del(OutputString)
    if dispPlots:
        import matplotlib.pyplot as plt
        
        
        width  = 6.5
        height = 5.5
        plt.figure(figsize=(width, height))
        
        plt.subplot(211)
        plt.grid()
        plt.ylabel('Angle of\n'+' Attack, deg', fontsize=20)
        # plt.ylim(0, 1.2*(180/pi)*max(alpha))
        plt.xlim(-0.5,0.5)
        plt.plot(y, (180/pi)*alphatw_vec,'k--', label=r'$\alpha_\mathrm{tw}$')
        plt.plot(y, (180/pi)*alpha_vec1,'k-', label=r'$\alpha_\mathrm{tw}+\alpha_w$')
        plt.yticks(fontsize=18)
        plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)),(), fontsize=18)
        plt.legend(fontsize=18,ncol = 2)
        
        plt.subplot(212)
        plt.grid()
        plt.ylabel(r'Vorticity, $\frac{\Gamma}{V_\infty b}$', fontsize=20)
        # plt.ylim(bottom=0)
        plt.xlim(-0.5,0.5)
        plt.plot(y, Gamma_vec1,'k-')
        # plt.plot(y, Gamma_vec2,'k-')
        plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)), fontsize=18)
        plt.yticks(fontsize=18)
        plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
        
        
        width  = 6.5
        height = 5.5
        plt.figure(figsize=(width, height))
        
        plt.subplot(212)
        plt.grid()
        plt.ylabel(r'Chord, $c/b$', fontsize=20)
        plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
        # plt.axis('equal')
        plt.plot(y, c_vec,'k-')
        plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)), fontsize=18)
        plt.yticks(fontsize=18)
        plt.ylim(0, 1.2*max(c_vec))
        plt.xlim(-0.5,0.5)
        
        width  = 6.5
        height = 5.5
        plt.figure(figsize=(width, height))
        
        plt.subplot(212)
        plt.grid()
        plt.xlabel(r'Wing Angle of Attack, deg', fontsize=20)
        plt.ylabel(r'Aspect Ratio', fontsize=20)
        dAR = AR_sweep[-1] - AR_sweep[0]
        yLims = [AR_sweep[-1] - dAR, AR_sweep[0]+dAR]
        LBoundArray = np.array([[alpha_w_Lims0[0], 1.00*yLims[0] + 0.00*yLims[1]],
                                [alpha_w_Lims0[0], 0.75*yLims[0] + 0.25*yLims[1]],
                                [alpha_w_Lims0[0], 0.50*yLims[0] + 0.50*yLims[1]],
                                [alpha_w_Lims0[0], 0.25*yLims[0] + 0.75*yLims[1]],
                                [alpha_w_Lims0[0], 0.00*yLims[0] + 1.00*yLims[1]]])
        UBoundArray = np.array([[alpha_w_Lims0[1], 1.00*yLims[0] + 0.00*yLims[1]],
                                [alpha_w_Lims0[1], 0.75*yLims[0] + 0.25*yLims[1]],
                                [alpha_w_Lims0[1], 0.50*yLims[0] + 0.50*yLims[1]],
                                [alpha_w_Lims0[1], 0.25*yLims[0] + 0.75*yLims[1]],
                                [alpha_w_Lims0[1], 0.00*yLims[0] + 1.00*yLims[1]]])
        plt.plot((180/pi)*alpha_w_sweep, AR_sweep,'k-')
        plt.plot((180/pi)*alpha_w, AR,'ko')
        plt.plot((180/pi)*LBoundArray[:,0], LBoundArray[:,1],'b>:')
        plt.plot((180/pi)*UBoundArray[:,0], UBoundArray[:,1],'r<:')
        plt.xticks(fontsize=18)
        plt.yticks(fontsize=18)
        plt.ylim(min(yLims), max(yLims))
        # plt.xlim((180/pi)*alpha_w_sweepLims[0], (180/pi)*alpha_w_sweepLims[-1])
        plt.xlim((180/pi)*min(alpha_w_sweepLims), (180/pi)*max(alpha_w_sweepLims))
        
        
        
        # # Plotting AVL Analysis
        # # creating grid for subplots
        # fig = plt.figure()
        # fig.set_figheight(3.5)
        # fig.set_figwidth(3.25)
         
        # ax1 = plt.subplot2grid(shape=(3, 1), loc=(0, 0), rowspan=2)
        # ax2 = plt.subplot2grid(shape=(3, 1), loc=(2, 0))
         
         
        # # initializing x,y axis value
        # x = np.arange(0, 10, 0.1)
        # y = np.cos(x)
         
        # # plotting subplots
        # ax1.plot(AVLdata1[:,1]/b, AVLdata1[:,4]/Cref,'tab:blue'  , linestyle='-' , label = r'$c_\ell c/c_\mathrm{ref}$')
        # ax1.plot(AVLdata1[:,1]/b, AVLdata1[:,7]     ,'tab:blue'  , linestyle='--', label = r'$c_\ell$')
        # ax1.plot(AVLdata2[:,1]/b, AVLdata2[:,4]/Cref,'tab:orange', linestyle='-' , label = r'$c_\ell c/c_\mathrm{ref}$')
        # ax1.plot(AVLdata2[:,1]/b, AVLdata2[:,7]     ,'tab:orange', linestyle='--', label = r'$c_\ell$')
        # ax1.set_xlim(-0.5,0.5)
        # ax1.set_ylim(0,1.2*max(AVLdata2[:,7]))
        # ax1.set_xticks([-0.5,-0.25,0,0.25,0.5])
        # ax1.grid(True)
        # ax1.legend(title=r'$c_L = $'+str(cL1)+'                 $c_L = $'+str(cL2), 
        #            fontsize=9, bbox_to_anchor=(0, 1, 1, 0), loc="lower left", mode="expand", ncol=2)
        # ax1.set_ylabel(r'$c_\ell$')
        

        # ax2.plot(AVLdata1[:,1]/b, -AVLdata1[:,5]*180/np.pi  ,'tab:blue'  , linestyle=':' , label = r'$c_\ell c/c_\mathrm{ref}$')
        # ax2.plot(AVLdata2[:,1]/b, -AVLdata2[:,5]*180/np.pi  ,'tab:orange', linestyle=':' , label = r'$c_\ell c/c_\mathrm{ref}$')
        # ax2.set_xlim(-0.5,0.5)
        # ax2.set_ylim(top=0)
        # ax2.set_ylim(-1.5*max(AVLdata2[:,5])*180/np.pi,0)
        # ax2.set_xticks([-0.5,-0.25,0,0.25,0.5])
        # ax2.grid(True)
        # # ax2.legend(fontsize=9, ncol = 2, loc = "lower center")
        # ax2.set_ylabel(r'$\alpha_i, deg$')
        # ax2.set_xlabel(r'$y/b$')
        
        # # automatically adjust padding horizontally
        # # as well as vertically.
        # plt.tight_layout()
         
        # # display plot
        # plt.show()
        
        
    
    # breakpoint()
    
    
    return [WingShape, l_vec2]