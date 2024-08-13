def TwoLift(n,b,AR,cL1,cL2,a,rho,W_G,liftShape1,liftShape2,dispPlots):
    
    """
    This function allows the user to specify two lift distributions
        distributions.
    
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
    a : Scalar [need to change to allow for distributions]
        airfoil lift curve slope
    rho : Scalar
        Air density
    W_G : Scalar
        Gross Aircraft Weight
    liftShape1 : 2 possible inputs
        'e'     can be given for an elliptical planform
    liftShape2 : nx2 numpy array
        'e'     can be given for an elliptical planform

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
    # b = b/1000
    
    S = b**2/AR
    
    # Define needed vectors and matrices ------------------------
    e = np.ones((n,1))
    I = np.eye(n)
    
    yF = np.linspace(0.5*(-1 + 1/n), 0.5*(1 - 1/n), n)
    y = np.zeros((n,1))
    for j in range(n):
        y[j,0] = yF[j]
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
    
    # Define Lift 1 Distribution ----------------------------------
    
    # eta = 2.68
    # # eta = 3.5
    if liftShape1[0] == 'eta' or liftShape2[0] == 'eta':
        K = np.zeros((n,n))
        for k in range(n):
            for j in range(n):
                K[j,k] = (1/n**2)*np.abs(j-k)
    
    if liftShape1 == 'e':
        l_sh1 = Qinv@e
    elif liftShape1 == 'p':
        l_sh1 = np.zeros((n,1))
        for j in range(n):
            l_sh1[j,0] = 1-4*y[j,0]**2
    elif liftShape1[0] == 'eta':
        eta1  = liftShape1[1]
        l_sh1 = Qinv@(I - eta1*K)@e
    else:
        print("!!!!!error in Lift Distribution 1!!!!!!!")
    
    if liftShape2 == 'e':
        l_sh2 = Qinv@e
    elif liftShape2[0] == 'eta':
        eta2  = liftShape2[1]
        l_sh2 = Qinv@(I - eta2*K)@e
    else:
        print("!!!!!error in Lift Distribution 2!!!!!!!")
    
    
    Gamma_vec1 = 0.5*(cL1/AR)*(n/(e.T@l_sh1))*l_sh1

    Gamma_vec2 = 0.5*(cL2/AR)*(n/(e.T@l_sh2))*l_sh2
    
    dGamma = Gamma_vec2 - Gamma_vec1
    
    
    # Iteratively solve for dalpha
    # cLalpha_0 = a/(1+(a/(np.pi*AR)))
    # dalpha_0 = (cL2 - cL1)/cLalpha_0
    
    dalpha_0 = ((n*(cL2-cL1)/(e.T@A@e))*(1 + (e.T@A@e)/(AR*n*pi)))[0,0]
    dalpha_Lims0 = np.array([max(2*Ainv@W@dGamma)[0],0])
    dalpha_LimDif = 2*AR*max(Ainv@dGamma)
    dalpha_Lims0[1] = dalpha_Lims0[0] + dalpha_LimDif
    print("dalpha_min = "+str(dalpha_Lims0[0]*180/pi)+"deg")
    print("dalpha_0   = "+str(    dalpha_0   *180/pi)+"deg")
    print("dalpha_max = "+str(dalpha_Lims0[1]*180/pi)+"deg")
    sweepLims = np.array([dalpha_Lims0[0] - dalpha_LimDif,dalpha_Lims0[1] + dalpha_LimDif])
    nSweep = 1001
    dAlphaSw = np.linspace(sweepLims[0],sweepLims[1],nSweep)
    ARsw = np.zeros(np.shape(dAlphaSw))
    for i in range(nSweep):
        ARsw[i] = n/(e.T@np.diagflat(1/(0.5*dAlphaSw[i]*A@e - W@dGamma))@dGamma)
        
    ARtol = 1e-7
    notFound=1
    count = 0
    dalpha_Lims = copy.deepcopy(dalpha_Lims0)
    while notFound and count < 100:
        count += 1
    #     print(count)
        
        dalpha = 0.5*(dalpha_Lims[0] + dalpha_Lims[1])
        
        M = np.diag(1/np.ndarray.flatten(0.5*A@(dalpha*e) - W@dGamma))
        ARcalc = (n/(e.T@M@dGamma))[0,0]
        
        ARres  = ((ARcalc - AR)/AR)
        # outputString = """$2.0i   dalpha = %6.4f;  AR = %6.4f;  AR_res = %6.4e;  """  % (count, dalpha*180/np.pi, ARcalc, ARres)
        print("""%2.0i   dalpha = %6.4f;  AR = %6.4f;  AR_res = %6.4e;  """  % (count, dalpha*180/np.pi, ARcalc, ARres))
        
        if abs(ARres) < ARtol:
            notFound = 0
        else:
            if ARres < 0:
                dalpha_Lims[0] = dalpha
            else:
                dalpha_Lims[1] = dalpha
    
    print("dalpha found!")
    
    # calculate chord and twist
    c_vec = (1/(0.5*dalpha*A@e - W@dGamma))*dGamma
    Cinv = np.diag(np.ndarray.flatten(1/c_vec))
    alpha_vec = 2*Ainv@(Cinv + W)@Gamma_vec1
    alpha_w = alpha_vec[n//2]
    alphatw_vec = alpha_vec - alpha_w
    alpha_vec2 = alpha_vec + dalpha
    
    # # Initial chord guess
    # # c0vec = (1/AR)*(n/(e.T@Qinv@e)[0,0])*Qinv@e
    # c0vec = (1/AR)*e
    
    # # Fixed-Point Iteration
    # maxIter = 100
    # import copy
    # cvec = np.zeros((n,maxIter))
    # cvec[:,0] = copy.deepcopy(c0vec[:,0])
    # C    = np.diag(np.ndarray.flatten(cvec[:,0]))
    # Cinv = np.diag(np.ndarray.flatten(1/cvec[:,0]))
    # dalpha = np.zeros(maxIter)

    # dalpha[0] = (1/(a*n))*(e.T@(Cinv + W)@(C@clvecD - 2*Gamma_vecD))[0,0];
    
    # eps_b = np.zeros(maxIter)
    # eps_c = np.zeros(maxIter)
    # converged = 0
    # eps_conv = 1e-10
    # i = 0
    # while converged == 0:
    #     i += 1

    #     c_u = Clinv@(a*dalpha[i-1]*inv(Cinv+W)@e + 2*Gamma_vecD)
    #     cvec[:,i] = omega*(n/AR)*(c_u/(e.T@c_u)[0,0])[:,0] + (1 - omega)*cvec[:,i-1]
        
    #     C    = np.diag(np.ndarray.flatten(cvec[:,i]))
    #     Cinv = np.diag(np.ndarray.flatten(1/cvec[:,i]))
    #     dalpha[i] = (1/(a*n))*(e.T@(Cinv + W)@(C@clvecD - 2*Gamma_vecD))[0,0]
        
    #     eps_b[i] = abs(dalpha[i] - dalpha[i - 1])
    #     eps_c[i] = (1/n)*e.T@abs(cvec[:,i] - cvec[:,i-1])
    #     outputString = "Aero Iter = %i ; eps_b = %4.3e  ; eps_c = %4.3e" \
    #         % (i,eps_b[i],eps_c[i])
    #     print(outputString)
    #     if eps_b[i] < eps_conv and eps_c[i] < eps_conv:
    #         converged = 1
    #     elif i == maxIter - 1:
    #         converged = 2
    
    # dalpha = dalpha[0:i+1]
    # cvec = cvec[:,0:i+1]
    # eps_b = eps_b[1:i]
    # eps_c = eps_c[1:i]
    
    # print("Fixed Point iteration complete!")
    # print("AR Check:",n/(e.T@cvec[:,i]))
    
    # alpha_vec1 = (2/a)*(Cinv + W)@Gamma_vecD
    # alphaw_1 = alpha_vec1[n//2]
    # alphaw_2 = alphaw_1 + dalpha[-1]
    # alphatw_vec = alpha_vec1 - alphaw_1
    # Gamma_vec1 = (a/2)*inv(Cinv + W)@(alphatw_vec+alphaw_1)
    # Gamma_vec2 = (a/2)*inv(Cinv + W)@(alphatw_vec+alphaw_2)
    # cl_vec1 = 2*Cinv@Gamma_vec1
    # cl_vec2 = 2*Cinv@Gamma_vec2
    # cL1 = 2*(AR/n)*(e.T@Gamma_vec1)[0,0]
    # cL2 = 2*(AR/n)*(e.T@Gamma_vec2)[0,0]
    
    # V1 = np.sqrt(2*W_G/(rho*S*cL1))
    # V2 = np.sqrt(2*W_G/(rho*S*cL2))
    
    # # l_vec1 = 0.5*rho*V1**2*C@cl_vec1
    # # l_vec2 = 0.5*rho*V2**2*C@cl_vec2
    
    # Lifting Line Math
#     cLalpha   = a*((e.T@inv(Cinv + W)@e)/(e.T@cvec[:,-1]))[0,0] #3D lift curve slope
#     cL0       = (a*AR/n)*(e.T@inv(Cinv+W)@alphatw_vec)[0,0]   #cL @ alpha_w = 0
    
#     OutputString ="""|------------------------------------------|
# |  cLalpha = %3.2f*a   ;    cL0 = %3.2f     |
# |------------------------------------------|
# | Condition  |   Lift Dist |    cl Dist    |
# | cL         |   %4.3f     |    %4.3f      |
# | V          |   %4.2f m/s |    %4.2f m/s  |
# | alpha_w    |   %4.3f deg |    %4.3f deg |
# |------------------------------------------|\n""" \
#     % (cLalpha/a,cL0,cL1, cL2, V1,V2,(180/pi)*alphaw_1 , (180/pi)*alphaw_2 )
#     print(OutputString)
#     del(OutputString)
    
    if dispPlots:
        import matplotlib.pyplot as plt
        
        width  = 6.5
        height = 6.5
        plt.figure(figsize=(width, height))
        plt.subplot(211)
        plt.grid()
        plt.ylabel(r'Vorticity, $\frac{\Gamma}{V_\infty b}$', fontsize=20)
        # plt.ylim(bottom=0)
        plt.xlim(-0.5,0.5)
        plt.plot(y, Gamma_vec1,'k--', label = r"$c_{L1}=$ %3.2f"%(cL1))
        plt.plot(y, Gamma_vec2,'k:', label = r"$c_{L2}=$ %3.2f"%(cL2))
        # plt.plot(y,     dGamma,'k:', linewidth=1, label = r"$\delta\Gamma$")
        plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)), fontsize=18)
        plt.yticks(fontsize=18)
        plt.ylim(top = 1.5*max(Gamma_vec2))
        plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
        plt.legend(fontsize=18,loc='upper center',ncol = 2)
        
        
        # " LL:cL=%3.2f; cDi=%5.4f"%(LLscalars[0,i],LLscalars[1,i])
        
        # plt.subplot(412)
        # plt.grid()
        # plt.ylabel('2D Lift \nCoefficient')
        # plt.plot(y, cl_vec1,'k-', linewidth =1, label = r'$@\alpha_1$')
        # plt.plot(y, cl_vec2,'k--', linewidth=1, label = r'$@\alpha_2$')
        # plt.plot(y, clvecD,'k:', linewidth=3, label = r"$c_{ l}*$")
        # # plt.plot(y, l_vec2*b/W_G,'k-', label = r'c_L = '+str(cL2))
        # plt.xlim(-0.5,0.5)
        # plt.xlabel(r'Spanwise Location, $y/b$', fontsize=10)
        # plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)), ())
        # plt.legend(fontsize=9,loc='center',ncol = 3)
        width  = 6.5
        height = 6.5
        plt.figure(figsize=(width, height))
        plt.subplot(211)
        plt.grid()
        plt.ylabel(r'Chord, $c/b$', fontsize=20)
        # plt.axis('equal')
        plt.plot(y, c_vec,'k-')
        plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)), ())
        plt.yticks(fontsize=18)
        # plt.yticks(np.array((0,0.05,0.10)))
        # plt.ylim(0, 1.2*max(cvec[:,-1]))
        plt.xlim(-0.5,0.5)
        
        plt.subplot(212)
        plt.grid()
        plt.ylabel(r'Angle of Attack, deg', fontsize=20)
        plt.ylim((180/pi)*min(alphatw_vec), 1.75*(180/pi)*max(alpha_vec2))
        plt.xlim(-0.5,0.5)
        # plt.plot(y, (180/pi)*alpha_vec2,'k:', label = r"$\vec{\alpha}$@$c_{L2}$")
        # plt.plot(y, (180/pi)*alpha_vec,'k--', label = r"$\vec{\alpha}$@$c_{L1}$")
        # plt.plot(y, (180/pi)*alphatw_vec,'k-', label = r"$\vec{\alpha}_{tw}$")
        plt.plot(y, (180/pi)*alpha_vec2,'k:', label = r"$@ c_{L2}$")
        plt.plot(y, (180/pi)*alpha_vec,'k--', label = r"$@ c_{L1}$")
        plt.plot(y, (180/pi)*alphatw_vec,'k-', label = r"$\vec{\alpha}_{tw}$")
        plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)),fontsize=16)
        plt.yticks(fontsize=18)
        # plt.yticks(np.array((0,-1,-2)))
        plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
        plt.legend(fontsize=18,loc='upper center',ncol = 3)
        
        width  = 6.5
        height = 6.5
        plt.figure(figsize=(width, height))
        plt.subplot(211)
        plt.grid()
        plt.ylabel(r'Aspect Ratio', fontsize=20)
        plt.xlabel(r'Diff. in Angle of Attack, deg', fontsize=20)
        # plt.axis('equal')
        plt.plot(dAlphaSw*180/pi, ARsw,'k-')
        plt.plot(dalpha*180/pi, ARcalc,'ko')
        # LBoundArray = np.array([[dalpha_Lims0[0],1.00*ARsw[-1] + 0.00*ARsw[0]],
        #                         [dalpha_Lims0[0],0.75*ARsw[-1] + 0.25*ARsw[0]],
        #                         [dalpha_Lims0[0],0.50*ARsw[-1] + 0.50*ARsw[0]],
        #                         [dalpha_Lims0[0],0.25*ARsw[-1] + 0.75*ARsw[0]],
        #                         [dalpha_Lims0[0],0.00*ARsw[-1] + 1.00*ARsw[0]]])
        # UBoundArray = np.array([[dalpha_Lims0[1],1.00*ARsw[-1] + 0.00*ARsw[0]],
        #                         [dalpha_Lims0[1],0.75*ARsw[-1] + 0.25*ARsw[0]],
        #                         [dalpha_Lims0[1],0.50*ARsw[-1] + 0.50*ARsw[0]],
        #                         [dalpha_Lims0[1],0.25*ARsw[-1] + 0.75*ARsw[0]],
        #                         [dalpha_Lims0[1],0.00*ARsw[-1] + 1.00*ARsw[0]]])
        # plt.plot((180/pi)*LBoundArray[:,0], LBoundArray[:,1],'b>:')
        # plt.plot((180/pi)*UBoundArray[:,0], UBoundArray[:,1],'r<:')
        plt.xticks(fontsize=18)
        plt.yticks(fontsize=18)
        plt.ylim(min(ARsw[0],ARsw[-1]), max(ARsw[0],ARsw[-1]))
        plt.xlim(dAlphaSw[0]*180/pi,dAlphaSw[-1]*180/pi)
        # plt.xlim(-0.5,0.5)
        
        # # breakpoint()
        # width  = 3.25
        # height = 2.5
        # plt.figure(figsize=(width, height))
        # plt.semilogy(np.linspace(1,i-1,i-1), eps_b,'k-', linewidth =1, label = r'$\epsilon_b$')
        # plt.semilogy(np.linspace(1,i-1,i-1), eps_c,'k--', linewidth =1, label = r'$\epsilon_c$')
        # plt.legend(fontsize=10,loc='lower center',ncol = 1)
        # # plt.xlim(-0.5,0.5)
        # plt.grid()
        # plt.ylabel(r'Convergence', fontsize=10)
        # plt.xlabel(r'Iteration', fontsize=10)
        # # plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)))     
    
    
    
    
    # c_vec = np.zeros((n,1))
    # for i in range(n):
    #     c_vec[i,0] = cvec[i,-1]
    # Return Wing Shape
    WingShape = np.concatenate((y,c_vec,alpha_vec),axis=1)    
    return WingShape


def TwoLiftCoefficient(n,b,AR,a,rho,W_G,clShape1,clShape2,dispPlots):
    
    """
    This function allows the user to specify two lift distributions
        distributions.
    
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
    a : Scalar [need to change to allow for distributions]
        airfoil lift curve slope
    rho : Scalar
        Air density
    W_G : Scalar
        Gross Aircraft Weight
    liftShape1 : 2 possible inputs
        'e'     can be given for an elliptical planform
    liftShape2 : nx2 numpy array
        'e'     can be given for an elliptical planform

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
    
    S = b**2/AR
    
    # Define needed vectors and matrices ------------------------
    e_vec = np.ones((n,1))
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
    Winv = ((2*pi)/(a*n))*Qinv
    
    A = a*I
    Ainv = (1/a)*I
    
    # Define Lift Coefficient Distributions -----------------------------------
    clShape1 = np.array(clShape1)
    clShape2 = np.array(clShape2)

    cl_vec1 = np.interp(y,clShape1[:,0],clShape1[:,1])
    cl_vec2 = np.interp(y,clShape2[:,0],clShape2[:,1])
    
    dcl_vec = cl_vec2 - cl_vec1
    Dclinv = np.diagflat(1/dcl_vec)
    
    dalpha = (((n/AR) + e_vec.T@Dclinv@Winv@dcl_vec)/(e_vec.T@Dclinv@Winv@A@e_vec))[0,0]
    c_vec = Dclinv@Winv@(dalpha*A@e_vec - dcl_vec)
    C = np.diagflat(c_vec)
    Cinv = np.diagflat(1/c_vec)
    
    cL1 = ((c_vec.T@cl_vec1)/(e_vec.T@c_vec))[0,0]
    cL2 = ((c_vec.T@cl_vec2)/(e_vec.T@c_vec))[0,0]
    
    alpha_vec1 = Ainv@(I + W@C)@cl_vec1
    alpha_w = alpha_vec1[n//2]
    alphatw_vec = alpha_vec1 - alpha_w
    alpha_vec2 = alpha_vec1 + dalpha
    
    Gamma_vec1 = (1/2)*inv(Cinv + W)@A@(alphatw_vec + alpha_w)
    
        
    if dispPlots:
        import matplotlib.pyplot as plt
        
        width  = 6.5
        height = 5.0
        plt.figure(figsize=(width, height))
        # plt.subplot(211)
        plt.grid()
        plt.ylabel(r'2D Lift Coeff.', fontsize=20)
        # plt.ylim(bottom=0)
        plt.xlim(-0.5,0.5)
        plt.plot(y, cl_vec1,'k--', label = r"$c_{\ell 1}$"%(cL1))
        plt.plot(y, cl_vec2,'k:', label = r"$c_{\ell 2}$"%(cL2))
        # plt.plot(y,     dGamma,'k:', linewidth=1, label = r"$\delta\Gamma$")
        plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)), fontsize=18)
        plt.yticks(fontsize=18)
        plt.ylim(top = 1.2*max(cl_vec2))
        plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
        plt.legend(fontsize=18,loc='upper center',ncol = 2)

        width  = 6.5
        height = 6.5
        plt.figure(figsize=(width, height))
        plt.subplot(211)
        plt.grid()
        plt.ylabel(r'Chord, $c/b$', fontsize=20)
        # plt.axis('equal')
        plt.plot(y, c_vec,'k-')
        plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)), ())
        plt.yticks(fontsize=18)
        # plt.yticks(np.array((0,0.05,0.10)))
        # plt.ylim(0, 1.2*max(cvec[:,-1]))
        plt.xlim(-0.5,0.5)
        
        plt.subplot(212)
        plt.grid()
        plt.ylabel(r'Angle of Attack, deg', fontsize=20)
        plt.ylim((180/pi)*min(alphatw_vec), 1.75*(180/pi)*max(alpha_vec2))
        plt.xlim(-0.5,0.5)
        # plt.plot(y, (180/pi)*alpha_vec2,'k:', label = r"$\vec{\alpha}$@$c_{L2}$")
        # plt.plot(y, (180/pi)*alpha_vec,'k--', label = r"$\vec{\alpha}$@$c_{L1}$")
        # plt.plot(y, (180/pi)*alphatw_vec,'k-', label = r"$\vec{\alpha}_{tw}$")
        plt.plot(y, (180/pi)*alpha_vec2,'k:', label = r"$@ \alpha_{w1}$")
        plt.plot(y, (180/pi)*alpha_vec1,'k--', label = r"$@ \alpha_{w1}$")
        plt.plot(y, (180/pi)*alphatw_vec,'k-', label = r"$\alpha_{tw}$")
        plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)),fontsize=16)
        plt.yticks(fontsize=18)
        # plt.yticks(np.array((0,-1,-2)))
        plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
        plt.legend(fontsize=18,loc='upper center',ncol = 3)
        
        
        # width  = 6.5
        # height = 5
        # plt.figure(figsize=(width, height))
        
        # ellshape = Qinv@e_vec
        # ellChord = (1/AR)*(n*ellshape/(e_vec.T@ellshape))
        
        # plt.subplot(211)
        # plt.grid()
        # plt.ylabel(r'$chord - ellipse$', fontsize=20)
        # plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
        # # plt.ylim(bottom=0)
        # plt.xlim(-0.5,0.5)
        # plt.plot(y, c_vec - ellChord,'k-')
        # plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)),fontsize=18)
        # plt.yticks(fontsize=20)
        
        # # breakpoint()
        # width  = 3.25
        # height = 2.5
        # plt.figure(figsize=(width, height))
        # plt.semilogy(np.linspace(1,i-1,i-1), eps_b,'k-', linewidth =1, label = r'$\epsilon_b$')
        # plt.semilogy(np.linspace(1,i-1,i-1), eps_c,'k--', linewidth =1, label = r'$\epsilon_c$')
        # plt.legend(fontsize=10,loc='lower center',ncol = 1)
        # # plt.xlim(-0.5,0.5)
        # plt.grid()
        # plt.ylabel(r'Convergence', fontsize=10)
        # plt.xlabel(r'Iteration', fontsize=10)
        # # plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)))     
    
    
    
    
    # c_vec = np.zeros((n,1))
    # for i in range(n):
    #     c_vec[i,0] = cvec[i,-1]
    # Return Wing Shape
    WingShape = np.concatenate((y,c_vec,alpha_vec1),axis=1)    
    return [WingShape,cL1,cL2]


def Lift_cl(n,b,AR,cLl,a,rho,W_G,liftShape,clShape,dispPlots):
    
    """
    This function allows the user to specify lift and lift coefficient
        distributions.
    
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
    
    S = b**2/AR
    
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
    
    A = a*I
    Ainv = (1/a)*I
    
    W = (n/(2*pi))*A@Q
    
    # Define Lift Distribution ----------------------------------
    if liftShape == 'e':
        l_sh = Qinv@e
        
    Vl         = np.sqrt(2*W_G/(rho*(b**2/AR)*cLl))   
    Gamma_vecD = (1/2)*(n*cLl/AR)*(l_sh/(e.T@l_sh))
    l_vecD     = rho*Vl*Gamma_vecD
    
    # Define cl Distribution
    clShape = np.array(clShape)
    clvecD = np.interp(y,clShape[:,0],clShape[:,1])
    
    Clinv = np.diag(np.ndarray.flatten(1/clvecD))
    clbar = (1/n)*(e.T@clvecD)[0,0]
    
    # Define relaxation
    if clbar < cLl:
        omega = clbar/cLl
    else:
        omega = 1
    
    # Initial chord guess
    # c0vec = (1/AR)*(n/(e.T@Qinv@e)[0,0])*Qinv@e
    c0vec = (1/AR)*e
    
    # Fixed-Point Iteration
    maxIter = 100
    import copy
    cvec = np.zeros((n,maxIter))
    cvec[:,0] = copy.deepcopy(c0vec[:,0])
    C    = np.diag(np.ndarray.flatten(cvec[:,0]))
    Cinv = np.diag(np.ndarray.flatten(1/cvec[:,0]))
    dalpha = np.zeros(maxIter)

    dalpha[0] = (1/(a*n))*(e.T@(Cinv + W)@(C@clvecD - 2*Gamma_vecD))[0,0];
    
    eps_d = np.zeros(maxIter)
    eps_c = np.zeros(maxIter)
    converged = 0
    eps_conv = 1e-10
    i = 0
    while converged == 0:
        i += 1

        c_u = Clinv@(a*dalpha[i-1]*inv(Cinv+W)@e + 2*Gamma_vecD)
        cvec[:,i] = omega*(n/AR)*(c_u/(e.T@c_u)[0,0])[:,0] + (1 - omega)*cvec[:,i-1]
        
        C    = np.diag(np.ndarray.flatten(cvec[:,i]))
        Cinv = np.diag(np.ndarray.flatten(1/cvec[:,i]))
        dalpha[i] = (1/(a*n))*(e.T@(Cinv + W)@(C@clvecD - 2*Gamma_vecD))[0,0]
        
        eps_d[i] = abs(dalpha[i] - dalpha[i - 1])
        eps_c[i] = (1/n)*e.T@abs(cvec[:,i] - cvec[:,i-1])
        outputString = "Aero Iter = %i ; eps_d = %4.3e  ; eps_c = %4.3e" \
            % (i,eps_d[i],eps_c[i])
        print(outputString)
        if eps_d[i] < eps_conv and eps_c[i] < eps_conv:
            converged = 1
        elif i == maxIter - 1:
            converged = 2
    
    dalpha = dalpha[0:i+1]
    cvec = cvec[:,0:i+1]
    eps_d = eps_d[1:i]
    eps_c = eps_c[1:i]
    
    print("Fixed Point iteration complete!")
    print("AR Check:",n/(e.T@cvec[:,i]))
    
    alpha_vec1 = (2/a)*(Cinv + W)@Gamma_vecD
    alphaw_1 = alpha_vec1[n//2]
    alphaw_2 = alphaw_1 + dalpha[-1]
    alphatw_vec = alpha_vec1 - alphaw_1
    Gamma_vec1 = (a/2)*inv(Cinv + W)@(alphatw_vec+alphaw_1)
    Gamma_vec2 = (a/2)*inv(Cinv + W)@(alphatw_vec+alphaw_2)
    cl_vec1 = 2*Cinv@Gamma_vec1
    cl_vec2 = 2*Cinv@Gamma_vec2
    cL1 = 2*(AR/n)*(e.T@Gamma_vec1)[0,0]
    cL2 = 2*(AR/n)*(e.T@Gamma_vec2)[0,0]
    
    V1 = np.sqrt(2*W_G/(rho*S*cL1))
    V2 = np.sqrt(2*W_G/(rho*S*cL2))
    
    e
    
    
    # l_vec1 = 0.5*rho*V1**2*C@cl_vec1
    # l_vec2 = 0.5*rho*V2**2*C@cl_vec2
    
    # Lifting Line Math
    cLalpha   = a*((e.T@inv(Cinv + W)@e)/(e.T@cvec[:,-1]))[0,0] #3D lift curve slope
    cL0       = (a*AR/n)*(e.T@inv(Cinv+W)@alphatw_vec)[0,0]   #cL @ alpha_w = 0
    
    OutputString ="""|------------------------------------------|
|  cLalpha = %3.2f*a   ;    cL0 = %3.2f     |
|------------------------------------------|
| Condition  |   Lift Dist |    cl Dist    |
| cL         |   %4.3f     |    %4.3f      |
| V          |   %4.2f m/s |    %4.2f m/s  |
| alpha_w    |   %4.3f deg |    %4.3f deg |
|------------------------------------------|\n""" \
    % (cLalpha/a,cL0,cL1, cL2, V1,V2,(180/pi)*alphaw_1 , (180/pi)*alphaw_2 )
    print(OutputString)
    del(OutputString)
    
    if dispPlots:
        import matplotlib.pyplot as plt
        
        width  = 6.5
        height = 6.5
        plt.figure(figsize=(width, height))
        
        plt.subplot(211)
        plt.grid()
        plt.ylabel(r'Vorticity, $\frac{\Gamma}{V_\infty b}$', fontsize=20)
        # plt.ylim(bottom=0)
        plt.xlim(-0.5,0.5)
        plt.ylim(top=1.5*max(max(Gamma_vec1),max(Gamma_vec2)))
        plt.plot(y, Gamma_vec1,'k--', linewidth = 3, label = r'$@\alpha_{w1}$ (Input)')
        plt.plot(y, Gamma_vec2, 'k:', linewidth = 1, label = r'$@\alpha_{w2}$')
        # plt.plot(y, Gamma_vecD,'k:', linewidth=3, label = r"$\vec{\Gamma}*$")
        plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)),())
        plt.yticks(fontsize=18)
        plt.legend(fontsize=18,loc='upper center',ncol = 2)
        
        plt.subplot(212)
        plt.grid()
        plt.ylabel('2D Lift \nCoefficient', fontsize = 20)
        plt.plot(y, cl_vec1,'k--', linewidth = 1, label = r'$@\alpha_{w1}$')
        plt.plot(y, cl_vec2, 'k:', linewidth = 3, label = r'$@\alpha_{w2}$ (Input)')
        # plt.plot(y, clvecD,'k:', linewidth=3, label = r"$\vec{c}_{\ell}*$")
        # plt.plot(y, l_vec2*b/W_G,'k-', label = r'c_L = '+str(cL2))
        plt.xlim(-0.5,0.5)
        # plt.xlabel(r'Spanwise Location, $y/b$', fontsize=10)
        plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)), fontsize = 18)
        plt.yticks(fontsize=18)
        plt.legend(fontsize=18,loc='center',ncol = 2)
        plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
        
        plt.figure(figsize=(width, height))
        plt.subplot(211)
        plt.grid()
        plt.ylabel(r'Chord, $c/b$', fontsize=20)
        # plt.axis('equal')
        # plt.plot(y, cvec[:,0],'k:')
        plt.plot(y, cvec[:,i-1],'k-')
        plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)), ())
        plt.yticks(fontsize=18)
        # plt.yticks(np.array((0,0.05,0.10)))
        # plt.ylim(0, 1.2*max(cvec[:,-1]))
        plt.xlim(-0.5,0.5)
        
        plt.subplot(212)
        plt.grid()
        plt.ylabel('Angle of\n'+'Attack (deg)', fontsize=20)
        # plt.ylim(0, 1.2*(180/pi)*max(alpha))
        plt.xlim(-0.5,0.5)
        plt.plot(y, (180/pi)*alphatw_vec           ,'k-', label = r'$\alpha_\mathrm{tw}$')
        plt.plot(y, (180/pi)*alphatw_vec + (180/pi)*alphaw_1,'k--',label = r'$@\alpha_{w1}$')
        plt.plot(y, (180/pi)*alphatw_vec + (180/pi)* alphaw_2,'k:', label = r'$@\alpha_{w2}$')
        plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)),fontsize=18)
        plt.yticks(fontsize=18)
        plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
        plt.legend(fontsize=18,loc='center', bbox_to_anchor=(0.5, 0.65),ncol = 3)
        
        # breakpoint()
        width  = 6.5
        height = 2.5
        plt.figure(figsize=(width, height))
        plt.semilogy(np.linspace(1,i-1,i-1), eps_d,'k--', label = r'$\epsilon_\delta$')
        plt.semilogy(np.linspace(1,i-1,i-1), eps_c,'k:', label = r'$\epsilon_c$')
        plt.legend(fontsize=18,loc='lower left',ncol = 1)
        # plt.xlim(-0.5,0.5)
        plt.grid()
        plt.ylabel('Convergence \n'+'Measure', fontsize=20)
        plt.xlabel(r'Iteration', fontsize=20)
        # iterPltMn  = 5
        # iterPltStp = 5
        # xRange     = np.ceil(i)/iterPltStp
        # iterTicks = np.linspace(iterPltMn, iterPltMn + xRange,iterPltStp)
        plt.xticks([5,10,15,20],fontsize=18)
        plt.yticks(fontsize=18)
        
        
    
    
    
    c_vec = np.zeros((n,1))
    for i in range(n):
        c_vec[i,0] = cvec[i,-1]
    # Return Wing Shape
    WingShape = np.concatenate((y,c_vec,alpha_vec1),axis=1)    
    return WingShape