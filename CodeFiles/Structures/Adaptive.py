def SimpleSparsWAil(WingShape,mb_vec,AFcoords,b,xcLL,sigma_max,FS,rho_m,g,w_e,dxLim,dxDfAM,sectBreaks,dispPlots,skinStructConsider = 1,UScurveOnly = 0):
    
    #import a whole buncha stuff
    import numpy as np
    from CodeFiles.Aerodynamics.Airfoils.operations import AFGeomAnalysis
    from CodeFiles.Aerodynamics.Airfoils.operations import AFBendProperties
    from CodeFiles.Aerodynamics.Airfoils.operations import AFCurveAnalysis
    from CodeFiles.Structures.functions import sectionSparPlace
    from CodeFiles.Structures.functions import wingSparAbs2NormPlace
    from CodeFiles.Structures.functions import wingSparAbsPlace
    from CodeFiles.Structures.functions import wingSparPlace
    from CodeFiles.Structures.functions import rSpars
    from CodeFiles.Structures.functions import bendingInertiaWAil
    from CodeFiles.Structures.functions import sectionBendingInertiaWAil
    from CodeFiles.Structures.functions import bisection    
    # breakpoint()
    
       
    # Define Thickness function of AF
    AFGeom = AFGeomAnalysis(AFcoords,101)
    tauDist = np.zeros((101,2))
    tauDist[:,0] = AFGeom[:,0]
    tauDist[:,1] = AFGeom[:,3]
    taumax  = max(tauDist[:,1])
    xtmax   = tauDist[np.argmax(tauDist[:,1]),0]
    
    # Initialize Spanwise locations and chord distribution
    n = len(WingShape[:,0])
    y_vec     = WingShape[:,0]
    c_vec = np.zeros((n,1))
    C     = np.zeros((n,n))
    for j in range(n):
        c_vec[j,0]  = b*WingShape[j,1] # NOTE c is now DIMENSIONAL
        C[j,j]      = b*WingShape[j,1]
    
    
    dxLim = dxLim/1000 # convert max spar spacing from mm to m
    
    ## ------------------------------------------------------------------- ##
    ## --------- Design Structure Based on Bending Inertia --------------- ##
    ## ------------------------------------------------------------------- ##
    # skinStructConsider = 0
    SkinStrength = 1.0
    sigma = sigma_max/FS
    
    [IxAFnorm,JzAFnorm,zcbar] = AFBendProperties(AFcoords,xcLL,dispPlots)
    IxSkin = np.zeros((n,1))
    for j in range(n):
        IxSkin[j,0] = SkinStrength*w_e*(c_vec[j,0])**3*IxAFnorm #skin bending intertia
        

    # IxReq_vec = Ireq_f_mb(c_vec,mb_vec,taumax,sigma)  
    IxReq_vec = (1/(2*sigma))*taumax*C@mb_vec
    IxReq_max = np.max(IxReq_vec)
    
    # Required functions for optimizer
    # costFunc    = lambda m_vec: (1/n)*sum(rSpars(c_vec,wingSparNormPlace(m_vec,xtmax),tauDist))[0]
    
    # breakpoint()
    
    # IxSparsMatLam  = lambda m_vec: bendingInertiaWAil(c_vec,wingSparNormPlace(m_vec,xtmax,dxLim/c_vec),tauDist,w_e,xcLL)
    IxSparsMatLam  = lambda m_vec: bendingInertiaWAil(c_vec,wingSparPlace(m_vec,xtmax,c_vec,xcLL,dxLim)[0],tauDist,w_e,xcLL)
    IxSparsVecLam  = lambda m_vec: [sum(IxSparsMatLam(m_vec)[i]) for i in range(n)]
    
    import copy 
    # Initial Guess of Structure
    IxMainSpar = (1/12)*w_e*(taumax*c_vec)**3
    if (IxMainSpar>IxReq_vec).all():
        m_vec0 = np.zeros((n,1))
        print("SKIN + MAIN SPAR STRONG ENOUGH STRUCTURE")
    else:
        iMaxStruct = np.argmax(IxReq_vec)
        # maxStructSolve = lambda m: 10**10*(sectionBendingInertiaWAil(c_vec[iMaxStruct,0],sectionSparPlace(m,m,xtmax,dxLim/c_vec[iMaxStruct,0]),tauDist,w_e,xcLL) - IxReq_max)
        if skinStructConsider == 1:
            maxStructSolve = lambda m: 10**10*(sectionBendingInertiaWAil(c_vec[iMaxStruct,0],sectionSparPlace(m,m,xtmax,dxLim/c_vec[iMaxStruct,0]),tauDist,w_e,xcLL) + IxSkin[iMaxStruct,0] - IxReq_max)
        else:
            maxStructSolve = lambda m: 10**10*(sectionBendingInertiaWAil(c_vec[iMaxStruct,0],sectionSparPlace(m,m,xtmax,dxLim/c_vec[iMaxStruct,0]),tauDist,w_e,xcLL) - IxReq_max)
        mMax = bisection(maxStructSolve,0,40,20)
        if mMax == None:
            print("STRUCTURE NOT POSSIBLE, LOAD TOO HIGH!!!!!!!")
        elif mMax == 0:
            m_vec0 = np.zeros((n,1))
        else:
            print('mMax = ',mMax)
            m_vec0 = np.zeros((n,1))
            # m_vec0 = copy.deepcopy(mb_vec)
            # if any(IxMainSpar > IxReq_vec):
            #     m_vec0[IxMainSpar > IxReq_vec] = max(m_vec0[IxMainSpar > IxReq_vec])
            # m_vec0 = m_vec0 - min(m_vec0)
            # # breakpoint()
            # m_vec0 = mMax*m_vec0/max(m_vec0)
    
    # m_vec       = copy.deepcopy(m_vec0)
    # xcSpars     = wingSparNormPlace(m_vec,xtmax)
    # xSpars   = wingSparAbsPlace(xcSpars,c_vec,xcLL)
    # IxSpars_vec = IxSparsVecLam(m_vec)
    
    # # Final Tweeks on Structure
    # for j in range(n):  #Sweep from left tip to right tip
    #     if IxReq_vec[j,0]>IxSpars_vec[j]:
    #         maxStructSolve = lambda m: 10**10*(sectionBendingInertiaWAil(c_vec[j,0],sectionSparPlace(m,xtmax),tauDist,w_e,xcLL) - IxReq_vec[j,0])
    #         m_vec[j,0] = bisection(maxStructSolve,0,20,20)
    m_vec = np.zeros((n,1))
    for j in range(n):  #Sweep from left tip to right tip
        if IxMainSpar[j,0]>IxReq_vec[j,0]:
            m_vec[j,0] = 0
        else:
            # structSolve = lambda m: 10**10*(sectionBendingInertiaWAil(c_vec[j,0],sectionSparPlace(m,xtmax,dxLim/c_vec[j,0]),tauDist,w_e,xcLL) - IxReq_vec[j,0])
            if skinStructConsider == 1:
                structSolve = lambda m: 10**10*(sectionBendingInertiaWAil(c_vec[j,0],sectionSparPlace(m,m,xtmax,dxLim/c_vec[j,0]),tauDist,w_e,xcLL) + IxSkin[j,0] - IxReq_vec[j,0])
            else:
                structSolve = lambda m: 10**10*(sectionBendingInertiaWAil(c_vec[j,0],sectionSparPlace(m,m,xtmax,dxLim/c_vec[j,0]),tauDist,w_e,xcLL) - IxReq_vec[j,0])
            m_vec[j,0] = bisection(structSolve,0,20,20)
            
            if  np.isnan(m_vec[j,0]):
                m_vec[j,0] = 0.0
        # print("j = ",j," ; m =", m_vec[j,0])
    
    # xcSpars      = wingSparNormPlace(m_vec,xtmax)
    # xSpars       = wingSparAbsPlace(xcSpars,c_vec,xcLL)
    print("Structure Solved")
    
    [xcSparsI,xSparsI,arc] = wingSparPlace(m_vec,xtmax,c_vec,xcLL,dxLim)
    print("Spars Placed")
    # IxSpars_vec  = IxSparsVecLam(m_vec)
    IxSpars_mat  = IxSparsMatLam(m_vec)
    IxSpars_vec  = [sum(IxSpars_mat[i]) for i in range(n)]
    print("Bending Inertia Calculated")
    # breakpoint()
    IxTot = IxSpars_vec + np.ndarray.flatten(IxSkin)
    
    Ixroot =  [IxTot[n//2], IxSkin[n//2,0], IxSpars_vec[n//2]]
    
    rSpars_vec   = rSpars(c_vec,xcSparsI,tauDist)
    wStruct_vec  = rho_m*g*w_e*rSpars_vec
    # breakpoint()
    FS_Spars_vec = FS*(IxSpars_vec/IxReq_vec[:,0])
    FS_Skin_vec  = FS*IxSkin[:,0]/IxReq_vec[:,0]
    # breakpoint()
    FS_Tot       = FS_Spars_vec + FS_Skin_vec
    FSroot = [FS_Tot[n//2], FS_Skin_vec[n//2], FS_Spars_vec[n//2],(IxSpars_vec+np.ndarray.flatten(IxSkin))[n//2] ,IxReq_max]
    print("Factor of Safety Calculated")
    
    
    ## --------------------------------------------------------------------- ##
    ## -------------- Limit Structure Base on Max Spacing ------------------ ##
    ## --------------------------------------------------------------------- ##
    
    print("Limiting Spar Spacing")
    nCurve = 50
    AFcurve = AFCurveAnalysis(AFcoords,nCurve)
    curveMapUS = np.zeros((nCurve,n))
    curveMapLS = np.zeros((nCurve,n))
    for j in range(n):
        curveMapUS[:,j] = AFcurve[:,1]/c_vec[j,0]
        curveMapLS[:,j] = AFcurve[:,2]/c_vec[j,0]
    
    nSurfPts = n*nCurve
    curvePts = np.zeros((nSurfPts,4))
    
    k = -1
    for j in range(n):
        for i in range(nCurve):
            k += 1
            curvePts[k,0] = y_vec[j]
            curvePts[k,1] = c_vec[j,0]*(AFcurve[i,0] - xcLL)
            curvePts[k,2] = curveMapUS[i,j]
            curvePts[k,3] = curveMapLS[i,j]
        
        
    ## Measure Spar Spacing for Inertially Placed Spars
    dxSparsI = np.zeros((n,np.shape(xSparsI)[1] - 1))
    for j in range(np.shape(xSparsI)[0]):
        for i in range((np.shape(xSparsI)[1] - 1)//2):
            if i == 0:
                dxSparsI[j,2*i    ] = -xSparsI[j,2*i + 1] + xSparsI[j,0]
                dxSparsI[j,2*i + 1] =  xSparsI[j,2*i + 2] - xSparsI[j,0]
            else:
                dxSparsI[j,2*i    ] = -xSparsI[j,2*i + 1] + xSparsI[j,2*i - 1]
                dxSparsI[j,2*i + 1] =  xSparsI[j,2*i + 2] - xSparsI[j,2*i]
    
    
    dxMaxWe   = dxDfAM[1]*w_e
    
    xSpars    = copy.deepcopy(xSparsI)
    xcSpars   = wingSparAbs2NormPlace(xSpars,c_vec,xcLL)
    kappaUS   = np.zeros(np.shape(xcSpars))
    kappaLS   = np.zeros(np.shape(xcSpars))
    
    
    dxSpars = np.zeros((n,np.shape(xSpars)[1] - 1))
    dxMax   = np.zeros(np.shape(dxSpars))
    dxMaxKappaUS   = np.zeros(np.shape(dxSpars))
    dxMaxKappaLS   = np.zeros(np.shape(dxSpars))
     
    
    
    for j in range(np.shape(xSpars)[0]):
        kappaUS[j,0] = np.interp(xcSpars[j,0], AFcurve[:,0], curveMapUS[:,j])
        kappaLS[j,0] = np.interp(xcSpars[j,0], AFcurve[:,0], curveMapLS[:,j])
    
    i = -1
    sparSpaceCheck = 1 
    while sparSpaceCheck == 1: # Sweeping through each spar
        i += 1
        print("Spar", i)
    
        for j in range(np.shape(xSpars)[0]): # Sweeping through spanwise loc
            if i == 0:
                try:
                    # breakpoint()
                    dxSpars[j,i] = -xSpars[j,i+1] + xSpars[j,i]
                except:
                    if j == 0:
                        xcSpars = np.append( xcSpars,np.zeros((n,1)),axis = 1)
                        xSpars  = np.append(  xSpars,wingSparAbsPlace(np.zeros((n,1)),c_vec,xcLL),axis = 1)
                        # breakpoint()
                        dxSpars = np.append( dxSpars,np.ones((n,1)),axis = 1)
                        dxMaxKappaUS = np.append( dxMaxKappaUS,np.zeros((n,1)),axis = 1)
                        dxMaxKappaLS = np.append( dxMaxKappaLS,np.zeros((n,1)),axis = 1)
                        dxMax        = np.append( dxMax       ,np.ones((n,1)),axis = 1)
                        # kappaUS = np.append( kappaUS,np.zeros((n,1)),axis = 1)
                        # kappaLS = np.append( kappaLS,np.zeros((n,1)),axis = 1)
                    # breakpoint  ()  
                    dxSpars[j,i] = -xSpars[j,i+1] + xSpars[j,i]
                # breakpoint()
                dxMaxKappaUS[j,i] = kappaUS[j,  0]*w_e**2/dxDfAM[0]
                dxMaxKappaLS[j,i] = kappaLS[j,  0]*w_e**2/dxDfAM[0]
                
                
            else:
                try:
                    # breakpoint()
                    if i%2 == 1: # aft spars
                        dxSpars[j,i] = xSpars[j, i + 1] - xSpars[j, i - 1]        
                    else: # forward spars
                        dxSpars[j,i] = xSpars[j, i - 1] - xSpars[j, i + 1]
                except:
                    if j == 0:
                        if i%2 == 1: # aft spars
                            xcSpars = np.append( xcSpars,np.ones((n,1)),axis = 1)
                            xSpars  = np.append(  xSpars,wingSparAbsPlace(np.ones((n,1)),c_vec,xcLL),axis = 1)
                        else: # forward spars
                            xcSpars = np.append( xcSpars,np.zeros((n,1)),axis = 1)   
                            xSpars  = np.append(  xSpars,wingSparAbsPlace(np.zeros((n,1)),c_vec,xcLL),axis = 1)
                        # breakpoint()
                        dxSpars = np.append( dxSpars,np.ones((n,1)),axis = 1)
                        dxMaxKappaUS = np.append( dxMaxKappaUS,np.zeros((n,1)),axis = 1)
                        dxMaxKappaLS = np.append( dxMaxKappaLS,np.zeros((n,1)),axis = 1)
                        dxMax        = np.append( dxMax       ,np.ones((n,1)),axis = 1)
                        # kappaUS = np.append( kappaUS,np.zeros((n,1)),axis = 1)
                        # kappaLS = np.append( kappaLS,np.zeros((n,1)),axis = 1)
                    # breakpoint  ()
                        if i%2 == 1: # aft spars
                            dxSpars[j,i] = xSpars[j, i + 1] - xSpars[j, i - 1]        
                        else: # forward spars
                            dxSpars[j,i] = xSpars[j, i - 1] - xSpars[j, i + 1]
                    
                dxMaxKappaUS[j,i] = kappaUS[j,i-1]*w_e**2/dxDfAM[0]
                dxMaxKappaLS[j,i] = kappaLS[j,i-1]*w_e**2/dxDfAM[0]
                    
            if UScurveOnly == 0:
                dxMax[j,i] = max(min(dxMaxKappaUS[j,i],dxMaxKappaLS[j,i]),dxMaxWe)
            else:
                dxMax[j,i] = max(dxMaxKappaUS[j,i],dxMaxWe)
            
            if dxSpars[j,i] < dxMax[j,i]:
                pass
            else:
                # print("Spar Edit")
                # breakpoint()
                dxSpars[j,i] = dxMax[j,i]
                if i == 0:
                    xSpars[j,i+1] = xSpars[j,i] - dxSpars[j,i]
                else:
                    if i%2 == 1: # aft spars
                        xSpars[j,i+1] = xSpars[j,i-1] + dxSpars[j,i]
                    else:
                        xSpars[j,i+1] = xSpars[j,i-1] - dxSpars[j,i]
            
            xcSpars[j,i+1] = wingSparAbs2NormPlace(xSpars[j,i+1],c_vec[j],xcLL)
            try:
                kappaUS[j,i+1] = np.interp(xcSpars[j,i+1], AFcurve[:,0], curveMapUS[:,j])
                kappaLS[j,i+1] = np.interp(xcSpars[j,i+1], AFcurve[:,0], curveMapLS[:,j])
            except:
                if j == 0:
                    kappaUS = np.append( kappaUS,np.zeros((n,1)),axis = 1)
                    kappaLS = np.append( kappaLS,np.zeros((n,1)),axis = 1)
                kappaUS[j,i+1] = np.interp(xcSpars[j,i+1], AFcurve[:,0], curveMapUS[:,j])
                kappaLS[j,i+1] = np.interp(xcSpars[j,i+1], AFcurve[:,0], curveMapLS[:,j])
        # breakpoint()
        if all(dxSpars[:,i] == 0) and all(dxSpars[:,i-1] == 0):
            sparSpaceCheck == 0
            # breakpoint()
            break
        if i > 200:
            print("INFITE LOOP DETECTED!")
            sparSpaceCheck == 0
            break
        
    
    
    # breakpoint()

    figReturn = 0
    if dispPlots:
        # AIAA Guidelines:
        #   3.25in wide, 3.5in tall
        #   min line weight: 0.28 pts
        #   min font size  : 8 pts
        
        import matplotlib.pyplot as plt
        print('I_x_skin_root  = %5.4e m^4' % IxSkin[n//2])
        print('I_x_spars_root = %5.4e m^4' % IxSpars_vec[n//2])
        print('I_x_root       = %5.4e m^4' % IxTot[n//2])
        print('frac_Iskin     = %5.4f m^4' % (IxSkin[n//2]/IxTot[n//2]))
        print('frac_Ispars    = %5.4f m^4' % (IxSpars_vec[n//2]/IxTot[n//2]))
        # m and bending inertia vs span
        width  = 6.5
        height = 7.0
        plt.figure(figsize=(width, height))
        
        plt.subplot(311)
        plt.grid()
        plt.ylabel('Spacing\n'+'Factor', fontsize=20)
        plt.xlim(-0.5,0.5)
        plt.plot(y_vec, m_vec,'k-',label = r'$m$')
        plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)), ())
        plt.yticks(fontsize=18)
        plt.ylim(None,1.2*max(m_vec))
        
        plt.subplot(212)
        plt.grid()
        plt.ylabel('Bending\n'+r'Inertia, mm$^4$', fontsize=20)
        plt.xlim(-0.5,0.5)
        # breakpoint()
        plt.plot(y_vec, [1000**4*IxSpars_vec[i] for i in range(n)],'k-',linewidth=1, label = 'Total')
        plt.plot(y_vec, 1000**4*IxReq_vec  ,'r--',linewidth=2, label = 'Req\'d')
        for i in range(len(IxSpars_mat[0,:])):
            if i == 0:
                # breakpoint()
                plt.plot(y_vec, 1000**4*IxSpars_mat[:,i],'b:',linewidth=1, label = 'Indiv.\n Spar')
            else:
                plt.plot(y_vec, 1000**4*IxSpars_mat[:,i],'b:',linewidth=1)
        
        # plt.ylim(0,1.5*1000**4*max(IxSpars_vec))
        plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)), fontsize=18)
        plt.yticks(fontsize=18)
        # plt.legend(fontsize=9,ncol=3,loc=3)
        plt.legend(fontsize=18,bbox_to_anchor=(0.02, 1.25, 0.95, .10), loc='upper center',
           ncol=3, mode="expand", borderaxespad=0.)
        plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
        
        # Factor of safety plots
        IxReq_vecPlot = np.array(IxReq_vec)*1000**4
        IxSpars_vecPlot = np.array(IxSpars_vec)*1000**4
        IxSkinPlot = np.array(IxSkin)*1000**4
        
        # Double plot
        width  = 6.5
        height = 6.5
        plt.figure(figsize=(width, height))
        plt.subplot(211)
        plt.grid()
        plt.xlim(-0.5,0.5)
        plt.plot(y_vec, IxReq_vecPlot  ,'-',color='r',linewidth=2, label = 'Req.')
        plt.plot(y_vec, IxSkinPlot     ,'-.',color=(0.9290, 0.6940, 0.1250),linewidth=2,label = r'Skin')
        plt.plot(y_vec, IxSpars_vecPlot,'--' ,color='b',linewidth=2,label = r'Spars')
        plt.ylim(None,1.5*max(IxSkinPlot))
        plt.yticks(fontsize=18)
        plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)), ())
        plt.ylabel('Bending\n'+'Inertia, $\mathrm{mm}^4$', fontsize=20)
        # plt.legend(fontsize=9,ncol=3,loc=3)
        plt.legend(fontsize=18,bbox_to_anchor=(0., 0.75, 1., .10), loc=3,
           ncol=4, mode="expand", borderaxespad=0.)
       
        plt.subplot(212)
        plt.grid()
        
        plt.xlim(-0.5,0.5)
        plt.plot(y_vec, FS_Tot      ,'k-',linewidth=2,label = r'Tot.')
        plt.plot([-0.5,0.5], [FS,FS],'r-',linewidth=2,label = r'Des.')
        plt.plot(y_vec, FS_Skin_vec ,'-.',color=(0.9290, 0.6940, 0.1250),linewidth=2,label = r'Skin')
        plt.plot(y_vec, FS_Spars_vec,'b--',linewidth=2,label = r'Spars')
        plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)),fontsize=18)
        plt.ylabel('Factor\n'+'of Safety', fontsize=20)
        plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
        plt.ylim(0,60)
        plt.yticks(fontsize=18)
        plt.legend(fontsize=18,bbox_to_anchor=(0.08, 0.9, .84, .10), loc=2,
           ncol=2, mode="expand", borderaxespad=0.)
        
        # Separate Plots
        width  = 6.5
        height = 3.5
        plt.figure(figsize=(width, height))
        plt.grid()
        plt.xlim(-0.5,0.5)
        plt.plot(y_vec, IxReq_vecPlot  ,'-',color='r',linewidth=2, label = 'Req.')
        plt.plot(y_vec, IxSkinPlot     ,'-.',color=(0.9290, 0.6940, 0.1250),linewidth=2,label = r'Skin')
        plt.plot(y_vec, IxSpars_vecPlot,'--' ,color='b',linewidth=2,label = r'Spars')
        plt.ylim(None,1.5*max(IxSkinPlot))
        plt.yticks(fontsize=18)
        plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)),fontsize=18)
        plt.ylabel('Bending\n'+'Inertia, $\mathrm{mm}^4$', fontsize=20)
        plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
        # plt.legend(fontsize=9,ncol=3,loc=3)
        plt.legend(fontsize=18,bbox_to_anchor=(0., 0.75, 1., .10), loc=3,
           ncol=4, mode="expand", borderaxespad=0.)
        
        # Plot Airfoil Only
        plt.figure(figsize=(width, height))
        plt.grid()
        
        plt.xlim(-0.5,0.5)
        plt.plot(y_vec, FS_Tot      ,'k-',linewidth=2,label = r'Tot.')
        plt.plot([-0.5,0.5], [FS,FS],'r-',linewidth=2,label = r'Des.')
        plt.plot(y_vec, FS_Skin_vec ,'-.',color=(0.9290, 0.6940, 0.1250),linewidth=2,label = r'Skin')
        plt.plot(y_vec, FS_Spars_vec,'b--',linewidth=2,label = r'Spars')
        plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)),fontsize=18)
        plt.ylabel('Factor\n'+'of Safety', fontsize=20)
        plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
        plt.ylim(0,60)
        plt.yticks(fontsize=18)
        plt.legend(fontsize=18,bbox_to_anchor=(0.08, 0.9, .84, .10), loc=2,
           ncol=2, mode="expand", borderaxespad=0.)
        
        # 
        width  = 6.5
        height = 2
        plt.figure(figsize=(width, height))
        
        plt.grid()
        plt.ylabel(r'$z/c$', fontsize=20)
        plt.xlabel(r'$x/c$', fontsize=20)
        plt.axis('equal')
        plt.plot(AFGeom[:,0], AFGeom[:,1], 'k-')
        plt.plot(AFGeom[:,0], AFGeom[:,2], 'k-')
        plt.plot(AFGeom[:,0], AFGeom[:,4], 'k--')
        plt.xticks(np.array((0,0.25,0.5,0.75,1.0)),fontsize = 18)
        plt.yticks(np.array((-0.2,-0.1,0.0,0.1,0.2)),fontsize = 18)
        # plt.ylim(0, 1.2*max(c_vec))
        plt.xlim(0.0,1.0)
        plt.ylim(-0.2,0.2)
        
        # Wing Top View Right Wing Only
        chordwiseShift = (-xcLL *c_vec/b)[n//2,0]
        width  = 6.5
        height = 3.5
        plt.figure(figsize=(width, height))
        # plt.grid()
        # for j in range(n):
        for i in range(len(xSparsI[0,:])):
            plt.plot(y_vec,  xSparsI[:,i]/b-chordwiseShift,'b-',linewidth=1)
        for i in range(len(sectBreaks)):
            cRib = np.interp(sectBreaks[i],y_vec,np.ndarray.flatten(c_vec))
            plt.plot([sectBreaks[i],sectBreaks[i]],  [-xcLL *cRib/b,(1-xcLL)*cRib/b]-chordwiseShift,'k-',linewidth=1)
        
        plt.plot([-0.5,0.5],  [-chordwiseShift,-chordwiseShift],'k--',linewidth=1.5)
        TE =    -xcLL*c_vec/b - chordwiseShift
        LE = (1-xcLL)*c_vec/b - chordwiseShift
        plt.plot(y_vec, TE,'k-',linewidth=1.5)
        plt.plot(y_vec, LE,'k-',linewidth=1.5)
        plt.plot([y_vec[-1], y_vec[-1]],[TE[-1], LE[-1]],'k-',linewidth=2)
        plt.xlim(0,0.5)
        plt.gca().invert_yaxis()
        plt.gca().set_aspect(1/b)
        plt.xlim(0.0,0.5)
        # plt.xticks(np.array((-0.50, -0.25, 0.0, 0.25, 0.50)))
        plt.xticks(np.array((0.0, 0.1, 0.2, 0.3, 0.4, 0.5)),fontsize=18)
        # plt.yticks(np.array((0.0, 0.075, 0.15)),fontsize=18)
        plt.yticks(fontsize=18)
        plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
        plt.ylabel(r'$x/b$', fontsize=20)        
        
        # Wing Top View Figure Split
        chordwiseShift = (-xcLL *c_vec/b)[n//2,0]
        width  = 6.5
        height = 7
        plt.figure(figsize=(width, height))
        plt.subplot(211)
        # plt.grid()
        # for j in range(n):
        for i in range(len(xSparsI[0,:])):
            plt.plot(y_vec,  xSparsI[:,i]/b-chordwiseShift,'b-',linewidth=1)
        for i in range(len(sectBreaks)):
            cRib = np.interp(sectBreaks[i],y_vec,np.ndarray.flatten(c_vec))
            plt.plot([sectBreaks[i],sectBreaks[i]],  [-xcLL *cRib/b,(1-xcLL)*cRib/b]-chordwiseShift,'k-',linewidth=1)
        
        plt.plot([-0.5,0.5],  [-chordwiseShift,-chordwiseShift],'k--',linewidth=1.5)
        TE =    -xcLL*c_vec/b - chordwiseShift
        LE = (1-xcLL)*c_vec/b - chordwiseShift
        plt.plot(y_vec, TE,'k-',linewidth=1.5)
        plt.plot(y_vec, LE,'k-',linewidth=1.5)
        plt.plot([y_vec[-1], y_vec[-1]],[TE[-1], LE[-1]],'k-',linewidth=2)
        plt.xlim(0,0.5)
        plt.gca().invert_yaxis()
        plt.gca().set_aspect(1/b)
        plt.xlim(0.0,0.5)
        # plt.xticks(np.array((-0.50, -0.25, 0.0, 0.25, 0.50)))
        plt.xticks(np.array((0.0, 0.1, 0.2, 0.3, 0.4, 0.5)),fontsize=18)
        # plt.yticks(np.array((0.0, 0.075, 0.15)),fontsize=18)
        plt.yticks(fontsize=18)
        plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
        plt.ylabel(r'$x/b$', fontsize=20)
        
        plt.subplot(212)
        # plt.grid()
        for i in range(len(xSparsI[0,:])):
            plt.plot(y_vec,  xSparsI[:,i]/b-chordwiseShift,'b-',linewidth=1)
        for i in range(len(sectBreaks)):
            cRib = np.interp(sectBreaks[i],y_vec,np.ndarray.flatten(c_vec))
            plt.plot([sectBreaks[i],sectBreaks[i]],  [-xcLL *cRib/b,(1-xcLL)*cRib/b]-chordwiseShift,'k-',linewidth=1)
        plt.plot([-0.5,0.5],  [-chordwiseShift,-chordwiseShift],'k--',linewidth=1.5)
        plt.plot(y_vec,  (-xcLL *c_vec/b)-chordwiseShift,'k-',linewidth=1.5)
        plt.plot(y_vec,((1-xcLL)*c_vec/b)-chordwiseShift,'k-',linewidth=1.5)
        plt.xlim(-0.5,0.0)
        plt.gca().invert_yaxis()
        plt.gca().set_aspect(1/b)
        plt.xticks(np.array((-0.5, -0.4, -0.3, -0.2, -0.1, 0.0)),fontsize=18)
        plt.yticks(fontsize=18)
        plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
        plt.ylabel(r'$x/b$', fontsize=20)
        
        # Wing Top View, Spars Only
        print("xcLL        = ",xcLL)
        print("c_vec[n//2] = ",c_vec[n//2])
        print("Chordshift  = ",chordwiseShift)
        
        chordwiseShift = (-xcLL *c_vec/b)[n//2,0]
        width  = 12
        height = 6
        plt.figure(figsize=(width, height))
        # plt.grid()
        for i in range(len(xSparsI[0,:])):
            plt.plot(y_vec,  xSparsI[:,i]/b-chordwiseShift,'b-',linewidth=1)
        # for i in range(len(sectBreaks)):
        #     cRib = np.interp(sectBreaks[i],y_vec,np.ndarray.flatten(c_vec))
        #     plt.plot([sectBreaks[i],sectBreaks[i]],  [-xcLL *cRib,(1-xcLL)*cRib]-chordwiseShift,'k-',linewidth=1)
        # plt.plot([-0.5,0.5],  [-chordwiseShift,-chordwiseShift],'k-',linewidth=2)
        plt.plot(y_vec,  -xcLL *c_vec/b-chordwiseShift,'k-',linewidth=1.5)
        plt.plot(y_vec,(1-xcLL)*c_vec/b-chordwiseShift,'k-',linewidth=1.5)
        plt.xlim(0,0.5)
        plt.xticks(np.array((-0.50, -0.25, 0.0, 0.25, 0.50)))
        plt.gca().invert_yaxis()
        plt.gca().set_aspect(1/b)
        plt.xticks(np.array((-0.5, -0.25, 0.0, 0.25, 0.5)),fontsize=18)
        plt.yticks(fontsize=18)
        plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
        plt.ylabel(r'$x/b$', fontsize=20)
        
        # Wing Top View, single plot
        print("xcLL        = ",xcLL)
        print("c_vec[n//2] = ",c_vec[n//2])
        print("Chordshift  = ",chordwiseShift)
        
        chordwiseShift = (-xcLL *c_vec/b)[n//2,0]
        width  = 12
        height = 6
        plt.figure(figsize=(width, height))
        # plt.grid()
        for i in range(len(xSparsI[0,:])):
            plt.plot(y_vec,  xSparsI[:,i]/b-chordwiseShift,'b-',linewidth=1)
        for i in range(len(sectBreaks)):
            cRib = np.interp(sectBreaks[i],y_vec,np.ndarray.flatten(c_vec))
            plt.plot([sectBreaks[i],sectBreaks[i]],  [-xcLL *cRib/b,(1-xcLL)*cRib/b]-chordwiseShift,'k-',linewidth=1)
        plt.plot([-0.5,0.5],  [-chordwiseShift,-chordwiseShift],'k--',linewidth=1.5)
        plt.plot(y_vec,  -xcLL *c_vec/b-chordwiseShift,'k-',linewidth=1.5)
        plt.plot(y_vec,(1-xcLL)*c_vec/b-chordwiseShift,'k-',linewidth=1.5)
        plt.xlim(0,0.5)
        plt.xticks(np.array((-0.50, -0.25, 0.0, 0.25, 0.50)))
        plt.gca().invert_yaxis()
        plt.gca().set_aspect(1/b)
        plt.xticks(np.array((-0.5, -0.25, 0.0, 0.25, 0.5)),fontsize=18)
        # plt.yticks(np.array((0.0, 0.075, 0.15)),fontsize=18)
        plt.yticks(fontsize=18)
        plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
        plt.ylabel(r'$x/b$', fontsize=20)
        
        
        # # Spar Spacing
        # chordwiseShift = (-xcLL *c_vec/b)[n//2,0]
        # width  = 6.5
        # height = 8
        # plt.figure(figsize=(width, height))
        
        # plt.subplot(411)
        # # plt.grid()
        # plt.plot(y_vec,  xSparsI[:,0]/b-chordwiseShift,'k--',linewidth=1)
        # for i in range(1,len(xSparsI[0,:])):
        #     plt.plot(y_vec,  xSparsI[:,i]/b-chordwiseShift,'--',linewidth=1)
        # for i in range(len(sectBreaks)):
        #     cRib = np.interp(sectBreaks[i],y_vec,np.ndarray.flatten(c_vec))
        #     plt.plot([sectBreaks[i],sectBreaks[i]],  [-xcLL *cRib/b,(1-xcLL)*cRib/b]-chordwiseShift,'k-',linewidth=1)
        # plt.plot([-0.5,0.5],  [-chordwiseShift,-chordwiseShift],'k-',linewidth=2)
        # plt.plot(y_vec,  -xcLL *c_vec/b-chordwiseShift,'k-',linewidth=2)
        # plt.plot(y_vec,(1-xcLL)*c_vec/b-chordwiseShift,'k-',linewidth=2)
        # plt.xlim(0,0.5)
        # plt.xticks(np.array((-0.50, -0.25, 0.0, 0.25, 0.50)))
        # plt.gca().invert_yaxis()
        # plt.gca().set_aspect(1/b)
        # plt.xticks(np.array((-0.5, -0.25, 0.0, 0.25, 0.5)),('','','','',''),fontsize=18)
        # # plt.yticks(np.array((0.0, 0.075, 0.15)),fontsize=18)
        # plt.yticks(fontsize=18)
        # # plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
        # plt.ylabel(r'$x/b$', fontsize=20)
        
        # plt.subplot(412)
        # # plt.grid()
        # for i in range(len(dxSparsI[0,:])):
        #     plt.plot(y_vec,  dxSparsI[:,i],'-',linewidth=1)
        # plt.xlim(0,0.5)
        # plt.xticks(np.array((-0.5, -0.25, 0.0, 0.25, 0.5)),('','','','',''),fontsize=18)
        # # plt.yticks(np.array((0.0, 0.075, 0.15)),fontsize=18)
        # plt.yticks(fontsize=18)
        # # plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
        # plt.ylabel(r'$dx$, m', fontsize=20)
        
        # colorList = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple',
        #              'tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
        # nColor = len(colorList)
        
        # plt.subplot(413)
        # # plt.grid()
        # plt.plot(y_vec,  xSpars[:,0]/b-chordwiseShift,'k-',linewidth=1)
        # for i in range(1,len(xSpars[0,:])):
        #     plt.plot(y_vec,  xSpars[:,i]/b-chordwiseShift,'-', color=colorList[i-1-int(np.floor(i/nColor))],linewidth=1)
        # for i in range(len(sectBreaks)):
        #     cRib = np.interp(sectBreaks[i],y_vec,np.ndarray.flatten(c_vec))
        #     plt.plot([sectBreaks[i],sectBreaks[i]],  [-xcLL *cRib/b,(1-xcLL)*cRib/b]-chordwiseShift,'k-',linewidth=1)
        # plt.plot([-0.5,0.5],  [-chordwiseShift,-chordwiseShift],'k--',linewidth=2)
        # plt.plot(y_vec,  -xcLL *c_vec/b-chordwiseShift,'k-',linewidth=2)
        # plt.plot(y_vec,(1-xcLL)*c_vec/b-chordwiseShift,'k-',linewidth=2)
        # plt.xlim(0,0.5)
        # plt.xticks(np.array((-0.50, -0.25, 0.0, 0.25, 0.50)))
        # plt.gca().invert_yaxis()
        # plt.gca().set_aspect(1/b)
        # plt.xticks(np.array((-0.5, -0.25, 0.0, 0.25, 0.5)),('','','','',''),fontsize=18)
        # plt.yticks(fontsize=18)
        # # plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
        # plt.ylabel(r'$x/b$', fontsize=20)
        
        # plt.subplot(414)
        # # plt.grid()
        # for i in range(len(dxSpars[0,:])-2):
        #     plt.plot(y_vec,  dxSpars[:,i],'-', color=colorList[i-int(np.floor(i/nColor))],linewidth=1)
        # plt.xlim(-0.5,0.5)
        # # plt.xticks(np.array((-0.50, -0.25, 0.0, 0.25, 0.50)))
        # # plt.gca().invert_yaxis()
        # # plt.gca().set_aspect(1/b)
        # plt.xticks(np.array((-0.5, -0.25, 0.0, 0.25, 0.5)),fontsize=18)
        # # plt.yticks(np.array((0.0, 0.075, 0.15)),fontsize=18)
        # plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
        # plt.ylabel(r'$dx$, m', fontsize=20)
        
        # Manufacturing-Constrained Wing Structure
        width  = 12
        height = 6
        plt.figure(figsize=(width, height))
        # plt.grid()
        plt.plot(y_vec,  xSpars[:,0]/b-chordwiseShift,'b-',linewidth=1)
        for i in range(1,len(xSpars[0,:])):
            plt.plot(y_vec,  xSpars[:,i]/b-chordwiseShift,'b-',linewidth=1)
        for i in range(len(sectBreaks)):
            cRib = np.interp(sectBreaks[i],y_vec,np.ndarray.flatten(c_vec))
            plt.plot([sectBreaks[i],sectBreaks[i]],  [-xcLL *cRib/b,(1-xcLL)*cRib/b]-chordwiseShift,'k-',linewidth=1)
        plt.plot([-0.5,0.5],  [-chordwiseShift,-chordwiseShift],'k--',linewidth=2)
        plt.plot(y_vec,  -xcLL *c_vec/b-chordwiseShift,'k-',linewidth=2)
        plt.plot(y_vec,(1-xcLL)*c_vec/b-chordwiseShift,'k-',linewidth=2)
        plt.xlim(0,0.5)
        plt.xticks(np.array((-0.50, -0.25, 0.0, 0.25, 0.50)))
        plt.gca().invert_yaxis()
        plt.gca().set_aspect(1/b)
        plt.xticks(np.array((-0.5, -0.25, 0.0, 0.25, 0.5)),fontsize=18)
        # plt.yticks(np.array((0.0, 0.075, 0.15)),fontsize=18)
        plt.yticks(fontsize=18)
        plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
        plt.ylabel(r'$x/b$', fontsize=20)
        
        # Manufacturing-Constrained Wing Structure (left wing only)
        width  = 6
        height = 4
        plt.figure(figsize=(width, height))
        # plt.grid()
        plt.plot(y_vec,  xSpars[:,0]/b-chordwiseShift,'b-',linewidth=1)
        for i in range(1,len(xSpars[0,:])):
            plt.plot(y_vec,  xSpars[:,i]/b-chordwiseShift,'b-',linewidth=1)
        for i in range(len(sectBreaks)):
            cRib = np.interp(sectBreaks[i],y_vec,np.ndarray.flatten(c_vec))
            plt.plot([sectBreaks[i],sectBreaks[i]],  [-xcLL *cRib/b,(1-xcLL)*cRib/b]-chordwiseShift,'k-',linewidth=1)
        plt.plot([-0.5,0],  [-chordwiseShift,-chordwiseShift],'k--',linewidth=2)
        plt.plot(y_vec,  -xcLL *c_vec/b-chordwiseShift,'k-',linewidth=2)
        plt.plot(y_vec,(1-xcLL)*c_vec/b-chordwiseShift,'k-',linewidth=2)
        plt.xlim(-0.5,0)
        plt.gca().invert_yaxis()
        plt.gca().set_aspect(1/b)
        plt.xticks(np.array((-0.5, -0.25, 0.0)),fontsize=18)
        # plt.yticks(np.array((0.0, 0.075, 0.15)),fontsize=18)
        plt.yticks(fontsize=18)
        plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
        plt.ylabel(r'$x/b$', fontsize=20)
        
        # Manufacturing-Constrained Wing Structure (Dimensional, right wing only)
        # yPlotMax = 660
        # width  = 6
        # height = 4
        # chordwiseShift = (-xcLL *c_vec/b)[n//2,0]
        # plt.figure(figsize=(width, height))
        # # plt.grid()
        # plt.plot(1000*b*y_vec,  1000*b*(xSpars[:,0]/b-chordwiseShift),'b-',linewidth=1)
        # for i in range(1,len(xSpars[0,:])):
        #     plt.plot(1000*b*y_vec, 1000*b*(xSpars[:,i]/b-chordwiseShift),'b-',linewidth=1)
        # for i in range(len(sectBreaks)):
        #     cRib = np.interp(sectBreaks[i],y_vec,np.ndarray.flatten(c_vec))
        #     plt.plot([1000*b*sectBreaks[i],1000*b*sectBreaks[i]],  [1000*(-xcLL)*cRib,1000*(1-xcLL)*cRib]-1000*b*chordwiseShift,'k-',linewidth=1)
        # plt.plot([0,1000*b*0.5],  [1000*b*-chordwiseShift,1000*b*-chordwiseShift],'k--',linewidth=2)
        # plt.plot(1000*b*y_vec, 1000*b*(-xcLL*c_vec/b-chordwiseShift),'k-',linewidth=2)
        # plt.plot(1000*b*y_vec, 1000*b*((1-xcLL)*c_vec/b-chordwiseShift),'k-',linewidth=2)
        # plt.plot([1000*b*0.5,1000*b*0.5], [1000*b*(-xcLL*c_vec[-1]/b-chordwiseShift),1000*b*((1-xcLL)*c_vec[-1]/b-chordwiseShift)],'k-',linewidth=2)
        
        # plt.gca().invert_yaxis()
        # plt.gca().set_aspect(1/b)
        # plt.xticks(fontsize=18)
        # # Iterating over all the axes in the figure 
        # # and make the Spines Visibility as False 
        # for pos in ['right', 'top', 'bottom']: 
        #         plt.gca().spines[pos].set_visible(False) 
        # # plt.yticks(np.array((0.0, 0.075, 0.15)),fontsize=18)
        # plt.yticks(fontsize=18)
        # plt.xlabel(r'Spanwise Location, $y$', fontsize=20)
        # plt.ylabel(r'$x$', fontsize=20)
        # plt.xlim(0,yPlotMax)
        
        # Manufacturing-Constrained Wing Structure (Dimensional, full wing)
        yPlotMax = 660
        width  = 12
        height = 4
        chordwiseShift = (-xcLL *c_vec/b)[n//2,0]
        plt.figure(figsize=(width, height))
        # plt.grid()
        plt.plot(1000*b*y_vec,  1000*b*(xSpars[:,0]/b-chordwiseShift),'b-',linewidth=1)
        for i in range(1,len(xSpars[0,:])):
            plt.plot(1000*b*y_vec, 1000*b*(xSpars[:,i]/b-chordwiseShift),'b-',linewidth=1)
        for i in range(len(sectBreaks)):
            cRib = np.interp(sectBreaks[i],y_vec,np.ndarray.flatten(c_vec))
            plt.plot([1000*b*sectBreaks[i],1000*b*sectBreaks[i]],  [1000*(-xcLL)*cRib,1000*(1-xcLL)*cRib]-1000*b*chordwiseShift,'k-',linewidth=1)
        plt.plot([-1000*b*0.5,1000*b*0.5],  [1000*b*-chordwiseShift,1000*b*-chordwiseShift],'k--',linewidth=2)
        plt.plot(1000*b*y_vec, 1000*b*(-xcLL*c_vec/b-chordwiseShift),'k-',linewidth=2)
        plt.plot(1000*b*y_vec, 1000*b*((1-xcLL)*c_vec/b-chordwiseShift),'k-',linewidth=2)
        plt.plot([1000*b*0.5,1000*b*0.5], [1000*b*(-xcLL*c_vec[-1]/b-chordwiseShift),1000*b*((1-xcLL)*c_vec[-1]/b-chordwiseShift)],'k-',linewidth=2)
        plt.plot([-1000*b*0.5,-1000*b*0.5], [1000*b*(-xcLL*c_vec[0]/b-chordwiseShift),1000*b*((1-xcLL)*c_vec[0]/b-chordwiseShift)],'k-',linewidth=2)
        
        plt.gca().invert_yaxis()
        plt.gca().set_aspect(1/b)
        plt.xticks(fontsize=18)
        # Iterating over all the axes in the figure 
        # and make the Spines Visibility as False 
        for pos in ['right', 'top', 'bottom']: 
                plt.gca().spines[pos].set_visible(False) 
        # plt.yticks(np.array((0.0, 0.075, 0.15)),fontsize=18)
        plt.yticks(fontsize=18)
        plt.xlabel(r'Spanwise Location $y$, mm', fontsize=20)
        plt.ylabel(r'$x$, mm', fontsize=20)
        plt.xlim(-yPlotMax,yPlotMax)
        
        if 1:
            # AF curvature
            width  = 6.5
            height = 4
            plt.figure(figsize=(width, height))
            
            plt.subplot(211)
            # plt.grid()
            plt.plot(AFcurve[:,0], AFcurve[:,1],'r-',linewidth=2)
            plt.plot(AFcurve[:,0], AFcurve[:,2],'b-',linewidth=2)
            plt.xlim(0,1)
            plt.ylim(0,3)
            # plt.xticks(np.array((-0.50, -0.25, 0.0, 0.25, 0.50)))
            # plt.gca().invert_yaxis()
            # plt.gca().set_aspect(1/b)
            # plt.xticks(np.array((-0.5, -0.25, 0.0, 0.25, 0.5)),fontsize=18)
            # plt.yticks(np.array((0.0, 0.075, 0.15)),fontsize=18)
            # plt.xlabel(r'Chordwise Location, $x/c$', fontsize=20)
            plt.ylabel(r'$\kappa c$', fontsize=20)
            # plt.title(r'Airfoil Curvature', fontsize=20)
            
            plt.subplot(212)
            # plt.grid()
            plt.plot(AFGeom[:,0], AFGeom[:,1],'r-',linewidth=2)
            plt.plot(AFGeom[:,0], AFGeom[:,2],'b-',linewidth=2)
            plt.xlim(0,1)
            # plt.ylim(0,3)
            # plt.xticks(np.array((-0.50, -0.25, 0.0, 0.25, 0.50)))
            # plt.gca().invert_yaxis()
            plt.gca().set_aspect(1)
            # plt.xticks(np.array((-0.5, -0.25, 0.0, 0.25, 0.5)),fontsize=18)
            # plt.yticks(np.array((0.0, 0.075, 0.15)),fontsize=18)
            plt.xlabel(r'Chordwise Location, $x/c$', fontsize=20)
            plt.ylabel(r'$z/b$', fontsize=20)
            
            # Spar Inspection Full wing
            chordwiseShift = (-xcLL *c_vec/b)[n//2,0]
            width  = 6.5
            height = 7
            plt.figure(figsize=(width, height))
            sparInspect = [2]
            
            plt.subplot(411)
            # plt.grid()
            # plt.plot(y_vec,  xSparsI[:,0]-chordwiseShift,'r--',linewidth=1)
            plt.plot(y_vec,  xSpars[:,0]/b-chordwiseShift,'k--',linewidth=1, label='Main Spar')
            for i in [sparInspect]:
                try:
                    plt.plot(y_vec,  xSparsI[:,i]/b-chordwiseShift,'r--',linewidth=1)
                except:
                    pass
                plt.plot(y_vec,  xSpars[:,i]/b-chordwiseShift,'k-',linewidth=1, label='1st Aft Spar')
            # for i in range(len(sectBreaks)):
            #     cRib = np.interp(sectBreaks[i],y_vec,np.ndarray.flatten(c_vec))
            #     plt.plot([sectBreaks[i],sectBreaks[i]],  [-xcLL *cRib,(1-xcLL)*cRib]-chordwiseShift,'k-',linewidth=1)
            # plt.plot([-0.5,0.5],  [-chordwiseShift,-chordwiseShift],'k-',linewidth=2)
            plt.plot(y_vec,  -xcLL *c_vec/b-chordwiseShift,'k-',linewidth=2)
            plt.plot(y_vec,(1-xcLL)*c_vec/b-chordwiseShift,'k-',linewidth=2)
            plt.xlim(-0.5,0.5)
            # plt.ylim(top = 0.2)
            plt.xticks(np.array((-0.50, -0.25, 0.0, 0.25, 0.50)), ['','','','',''])
            plt.gca().invert_yaxis()
            plt.gca().set_aspect(1/b)
            # plt.xticks(np.array((-0.5, -0.25, 0.0, 0.25, 0.5)),fontsize=18)
            # plt.yticks(np.array((0.0, 0.075, 0.15)),fontsize=18)
            plt.yticks(fontsize=18)
            # plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
            plt.ylabel(r'$x/b$', fontsize=20)
            # plt.legend(loc='lower center',ncol = 2,fontsize=16)
            plt.legend(fontsize=18,bbox_to_anchor=(0.02, 1.25, 0.95, .10), loc='lower center',
               ncol=2, mode="expand", borderaxespad=0.)
            
            plt.subplot(412)
            # plt.grid()
            # for i in range(len(dxSparsI[0,:])):
            for i in sparInspect:
                plt.semilogy(     y_vec,  kappaUS[:,i-2], 'r--',linewidth=1,label=r"US")
                plt.semilogy(     y_vec,  kappaLS[:,i-2], 'b--',linewidth=1,label=r"LS")
            plt.xlim(-0.5,0.5)
            # plt.xticks(np.array((-0.50, -0.25, 0.0, 0.25, 0.50)))
            # plt.gca().invert_yaxis()
            # plt.gca().set_aspect(1/b)
            plt.xticks(np.array((-0.50, -0.25, 0.0, 0.25, 0.50)), ['','','','',''])
            plt.yticks(fontsize=18)
            # plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
            plt.ylabel(r'$\kappa$, 1/m', fontsize=20)
            plt.legend(loc='upper center',ncol = 2,fontsize=16)
            
            plt.subplot(212)
            # plt.grid()
            # for i in range(len(dxSparsI[0,:])):
            for i in sparInspect:
                # breakpoint()
                plt.plot(     y_vec,         1000*dxMax[:,i-1], 'C7--',linewidth=2.5,label=r"active")
                plt.plot(     y_vec,       1000*dxSpars[:,i-1], 'C7:' ,linewidth=2,label=r"built")
                plt.plot(     y_vec,  1000*dxMaxKappaUS[:,i-1], 'r-',linewidth=1,label=r"$\kappa,US$")
                plt.plot(     y_vec,  1000*dxMaxKappaLS[:,i-1], 'b-',linewidth=1,label=r"$\kappa,LS$")
                plt.plot([-0.5,0.5],  [1000*dxMaxWe,1000*dxMaxWe], 'g-',linewidth=1,label=r"$w_e$")
            plt.xlim(-0.5,0.5)
            plt.ylim( 0.0,2.0*1000*max(dxSpars[:,sparInspect[0]-1]))
            plt.xticks(np.array((-0.5, -0.25, 0.0, 0.25, 0.5)),fontsize=18)
            plt.yticks(fontsize=18)
            plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
            plt.ylabel(r'Spar Spacing, mm', fontsize=20)
            plt.legend(loc='upper center',ncol = 3,fontsize=16)
            # plt.legend(loc='best',ncol = 3,fontsize=16)
   
        # # Spar Inspection Right Wing ONLY
        # cb_vec = c_vec/b
        # chordwiseShift = (-xcLL *cb_vec)[n//2,0]
        # width  = 6.5
        # height = 10
        # plt.figure(figsize=(width, height))
        # sparInspect = [2]
        
        # plt.subplot(311)
        # # plt.grid()
        # plt.plot(y_vec,  xSparsI[:,0]-chordwiseShift,'r--',linewidth=1)
        # plt.plot(y_vec,  xSpars[:,0]-chordwiseShift,'k--',linewidth=1)
        # for i in [sparInspect]:
        #     try:
        #         plt.plot(y_vec,  xSparsI[:,i]-chordwiseShift,'r--',linewidth=1)
        #     except:
        #         pass
        #     plt.plot(y_vec,  xSpars[:,i]-chordwiseShift,'k--',linewidth=1)
        # # for i in range(len(sectBreaks)):
        # #     cRib = np.interp(sectBreaks[i],y_vec,np.ndarray.flatten(c_vec))
        # #     plt.plot([sectBreaks[i],sectBreaks[i]],  [-xcLL *cRib,(1-xcLL)*cRib]-chordwiseShift,'k-',linewidth=1)
        # plt.plot([-0.5,0.5],  [-chordwiseShift,-chordwiseShift],'k-',linewidth=2)
        # plt.plot(y_vec,  -xcLL *cb_vec-chordwiseShift,'k-',linewidth=2)
        # plt.plot(y_vec,(1-xcLL)*cb_vec-chordwiseShift,'k-',linewidth=2)
        # plt.xlim(0,0.5)
        # plt.gca().invert_yaxis()
        # plt.gca().set_aspect(1/b)
        # plt.xticks(np.array((0.0, 0.25, 0.50)), ['','',''])
        # plt.yticks(fontsize=18)
        # # plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
        # plt.ylabel(r'$x/b$', fontsize=20)
        
        # plt.subplot(312)
        # # plt.grid()
        # # for i in range(len(dxSparsI[0,:])):
        # for i in [sparInspect]:
        #     plt.semilogy(     y_vec,  kappaUS[:,i], 'r--',linewidth=1)
        #     plt.semilogy(     y_vec,  kappaLS[:,i], 'b--',linewidth=1)
        # plt.xlim(0,0.5)
        # # plt.xticks(np.array((-0.50, -0.25, 0.0, 0.25, 0.50)))
        # # plt.gca().invert_yaxis()
        # # plt.gca().set_aspect(1/b)
        # plt.xticks(np.array((0.0, 0.25, 0.50)), ['','',''])
        # plt.yticks(fontsize=18)
        # # plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
        # plt.ylabel(r'$\kappa$, 1/m', fontsize=20)
        
        # plt.subplot(313)
        # # plt.grid()
        # # for i in range(len(dxSparsI[0,:])):
        # for i in sparInspect:
        #     # breakpoint()
        #     plt.plot(     y_vec,  dxMaxKappaUS[:,i-1], 'r--',linewidth=1,label=r"$dx_{\kappa,US}$")
        #     plt.plot(     y_vec,  dxMaxKappaLS[:,i-1], 'b--',linewidth=1,label=r"$dx_{\kappa,LS}$")
        #     plt.plot([-0.5,0.5],  [dxMaxWe,dxMaxWe], 'g--',linewidth=1,label=r"$dx_{\kappa,We}$")
        #     plt.plot(     y_vec,         dxMax[:,i-1], 'k--',linewidth=1,label=r"$dx_{max}$")
        #     plt.plot(     y_vec,       dxSpars[:,i-1], 'k-' ,linewidth=1,label=r"$dx$")
        # plt.xlim(0.0,0.5)
        # # breakpoint()
        # plt.ylim( 0.0,1.5*max(max(dxSpars[:,sparInspect[0]-1]),min(dxMaxKappaUS[:,i-1])))
        # # plt.xticks(np.array((-0.50, -0.25, 0.0, 0.25, 0.50)))
        # # plt.gca().invert_yaxis()
        # # plt.gca().set_aspect(1/b)
        # plt.xticks(np.array((0.0, 0.25, 0.5)),fontsize=18)
        # plt.yticks(fontsize=18)
        # plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
        # plt.ylabel(r'$dx$, m', fontsize=20)
        # plt.legend(loc=2,ncol = 3, fontsize=17)
        
        
        
        # # Skin Analysis
        # width  = 6.5
        # height = 9
        # plt.figure(figsize=(width, height))
        
        # plt.subplot(411)
        # CS = plt.tricontourf(curvePts[:,0], curvePts[:,1]-chordwiseShift, curvePts[:,2],[0,1,2,3,5,10], cmap="YlOrRd_r")
        # plt.colorbar(CS, orientation = "horizontal")
        # # for i in range(len(xSparsAM[0,:])):
        # #     plt.plot(y_vec,  xSparsAM[:,i]-chordwiseShift,'k--',linewidth=1)
        # # for i in range(len(sectBreaks)):
        # #     cRib = np.interp(sectBreaks[i],y_vec,np.ndarray.flatten(c_vec))
        # #     plt.plot([sectBreaks[i],sectBreaks[i]],  [-xcLL *cRib,(1-xcLL)*cRib]-chordwiseShift,'k-',linewidth=1)
        # # plt.plot([-0.5,0.5],  [-chordwiseShift,-chordwiseShift],'k-',linewidth=1)
        # plt.plot(y_vec,  -xcLL *c_vec-chordwiseShift,'k-',linewidth=2)
        # plt.plot(y_vec,(1-xcLL)*c_vec-chordwiseShift,'k-',linewidth=2)
        # plt.xlim(0,0.5)
        # plt.xticks(np.array((-0.50, -0.25, 0.0, 0.25, 0.50)))
        # plt.gca().invert_yaxis()
        # plt.gca().set_aspect("equal")
        # plt.xlim(-0.5,0.5)
        # plt.xlabel(r'Spanwise Location, [y/b]', fontsize=10)
        # plt.ylabel('Chordwise \nLocation\n'+r'[x/b]')
        # plt.title('Upper Surface Curvature')
        
        
        # plt.subplot(412)
        # CS = plt.tricontourf(curvePts[:,0], curvePts[:,1]-chordwiseShift, curvePts[:,3],[0,1,2,3,5,10], cmap="YlOrRd_r")
        # plt.colorbar(CS, orientation = "horizontal")
        # # for i in range(len(xSparsAM[0,:])):
        # #     plt.plot(y_vec,  xSparsAM[:,i]-chordwiseShift,'k--',linewidth=1)
        # # for i in range(len(sectBreaks)):
        # #     cRib = np.interp(sectBreaks[i],y_vec,np.ndarray.flatten(c_vec))
        # #     plt.plot([sectBreaks[i],sectBreaks[i]],  [-xcLL *cRib,(1-xcLL)*cRib]-chordwiseShift,'k-',linewidth=1)
        # # plt.plot([-0.5,0.5],  [-chordwiseShift,-chordwiseShift],'k-',linewidth=1)
        # plt.plot(y_vec,  -xcLL *c_vec-chordwiseShift,'k-',linewidth=2)
        # plt.plot(y_vec,(1-xcLL)*c_vec-chordwiseShift,'k-',linewidth=2)
        # plt.xlim(0,0.5)
        # plt.xticks(np.array((-0.50, -0.25, 0.0, 0.25, 0.50)))
        # plt.gca().invert_yaxis()
        # plt.gca().set_aspect("equal")
        # plt.xlim(-0.5,0.5)
        # plt.xlabel(r'Spanwise Location, [y/b]', fontsize=10)
        # plt.ylabel('Chordwise \nLocation\n'+r'[x/b]')
        # plt.title('Lower Surface Curvature')
        
        # plt.subplot(413)
        # # CS = plt.tricontourf(curvePts[:,0], curvePts[:,1]-chordwiseShift, curvePts[:,3],[0,1,2,3,5,10], cmap="YlOrRd_r")
        # # plt.colorbar(CS, orientation = "horizontal")
        # for i in range(len(xSparsAM[0,:])):
        #     plt.plot(y_vec,  xSparsAM[:,i]-chordwiseShift,'k--',linewidth=1)
        # for i in range(len(sectBreaks)):
        #     cRib = np.interp(sectBreaks[i],y_vec,np.ndarray.flatten(c_vec))
        #     plt.plot([sectBreaks[i],sectBreaks[i]],  [-xcLL *cRib,(1-xcLL)*cRib]-chordwiseShift,'k-',linewidth=1)
        # plt.plot([-0.5,0.5],  [-chordwiseShift,-chordwiseShift],'k-',linewidth=1)
        # plt.plot(y_vec,  -xcLL *c_vec-chordwiseShift,'k-',linewidth=2)
        # plt.plot(y_vec,(1-xcLL)*c_vec-chordwiseShift,'k-',linewidth=2)
        # plt.xlim(0,0.5)
        # plt.xticks(np.array((-0.50, -0.25, 0.0, 0.25, 0.50)))
        # plt.gca().invert_yaxis()
        # plt.gca().set_aspect("equal")
        # plt.xlim(-0.5,0.5)
        # # plt.xlabel(r'Spanwise Location, [y/b]', fontsize=10)
        # plt.ylabel('Chordwise \nLocation\n'+r'[x/b]')
        # # plt.title('Lower Surface Curvature')
        
        # plt.subplot(414)
        # plt.grid()
        # plt.ylabel('Unsupported \n Arc Length [mm]', fontsize=10)
        # plt.xlim(-0.5,0.5)
        # plt.plot(y_vec, 1000*arc[:,0],'k:' ,label = r'$dx_\mathrm{fore}$')
        # plt.plot(y_vec, 1000*arc[:,1],'k--',label = r'$dx_\mathrm{aft}$')
        # plt.ylim(0,1500*max(max(arc[:,0]),max(arc[:,1])))
        # # plt.ylim(None,1.5*max(IxSpars_vec))
        # plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)))
        # # plt.legend(fontsize=9,ncol=3,loc=3)
        # plt.xlabel(r'Spanwise Location, [y/b]', fontsize=10)
        # plt.legend(fontsize=9,bbox_to_anchor=(0., 0.75, 1., .10), loc=3,
        #    ncol=4, mode="expand", borderaxespad=0.)
        # plt.title('Spar Spacing')
        
        
        # JustinNumber = np.zeros((nSurfPts))
        # k = -1
        # for j in range(n):
        #     for i in range(nCurve):
        #         k += 1
        #         # curvePts[k,0] = y_vec[j]
        #         # curvePts[k,1] = c_vec[j,0]*(AFcurve[i,0] - xcLL)
        #         # curvePts[k,2] = curveMapUS[i,j]
        #         if AFcurve[i,0] < xcLL:
        #             dx = 1000*arc[j,0]
        #         else:
        #             dx = 1000*arc[j,1]
        #         JustinNumber[k] = curvePts[k,2]*dx
        
        # # plt.subplot(313)
        # # CS = plt.tricontourf(curvePts[:,0], curvePts[:,1]-chordwiseShift, JustinNumber,np.linspace(10,200,8), cmap="coolwarm")
        # # plt.colorbar(CS,orientation = "horizontal")
        # # for j in range(n):
        # #     for i in range(len(xSpars[j])):
        # #         plt.plot(y_vec[j],  xSpars[j][i]-chordwiseShift,'k.',markersize=0.75)
        # # for i in range(len(sectBreaks)):
        # #     cRib = np.interp(sectBreaks[i],y_vec,np.ndarray.flatten(c_vec))
        # #     plt.plot([sectBreaks[i],sectBreaks[i]],  [-xcLL *cRib,(1-xcLL)*cRib]-chordwiseShift,'k-',linewidth=1)
        # # plt.plot([-0.5,0.5],  [-chordwiseShift,-chordwiseShift],'k-',linewidth=1)
        # # plt.plot(y_vec,  -xcLL *c_vec-chordwiseShift,'k-',linewidth=2)
        # # plt.plot(y_vec,(1-xcLL)*c_vec-chordwiseShift,'k-',linewidth=2)
        # # plt.xlim(0,0.5)
        # # plt.xticks(np.array((-0.50, -0.25, 0.0, 0.25, 0.50)))
        # # plt.gca().invert_yaxis()
        # # plt.gca().set_aspect("equal")
        # # plt.xlim(-0.5,0.5)
        # # plt.xlabel(r'Spanwise Location, [y/b]', fontsize=10)
        # # plt.ylabel('Chordwise \nLocation\n'+r'[x/b]')
        
        # plt.subplot(313)
        # CS = plt.tricontourf(curvePts[:,0], curvePts[:,1]-chordwiseShift, curvePts[:,3],[0,1,2,3,5,10], cmap="YlOrRd_r")
        # plt.colorbar(CS,orientation = "horizontal")
        # for j in range(n):
        #     for i in range(len(xSpars[j])):
        #         plt.plot(y_vec[j],  xSpars[j][i]-chordwiseShift,'k.',markersize=0.75)
        # for i in range(len(sectBreaks)):
        #     cRib = np.interp(sectBreaks[i],y_vec,np.ndarray.flatten(c_vec))
        #     plt.plot([sectBreaks[i],sectBreaks[i]],  [-xcLL *cRib,(1-xcLL)*cRib]-chordwiseShift,'k-',linewidth=1)
        # plt.plot([-0.5,0.5],  [-chordwiseShift,-chordwiseShift],'k-',linewidth=1)
        # plt.plot(y_vec,  -xcLL *c_vec-chordwiseShift,'k-',linewidth=2)
        # plt.plot(y_vec,(1-xcLL)*c_vec-chordwiseShift,'k-',linewidth=2)
        # plt.xlim(0,0.5)
        # plt.xticks(np.array((-0.50, -0.25, 0.0, 0.25, 0.50)))
        # plt.gca().invert_yaxis()
        # plt.gca().set_aspect("equal")
        # plt.xlim(-0.5,0.5)
        # plt.xlabel(r'Spanwise Location, [y/b]', fontsize=10)
        # plt.ylabel('Chordwise \nLocation\n'+r'[x/b]')
        # plt.title('Lower Surface Curvature')
        
        
        
        
        
        # # Test Section
        # testSection = 6
        # width  = 2.5
        # height = 4
        # plt.figure(figsize=(width, height))
        # breakpoint()
        # ySect = [b*sectBreaks[testSection-1],b*sectBreaks[testSection]]
        # yDist = b*y_vec - ySect[0]
        # yMid = 0.5*(ySect[0]+ySect[1])
        # plt.subplot(2,1,1)
        # plt.grid()
        # plt.ylabel('Unsupported \n Wall Length [mm]', fontsize=10)
        # plt.plot(1000*yDist, 1000*arc[:,0],'k:' ,label = r'$dx_\mathrm{fore}$')
        # plt.plot(1000*yDist, 1000*arc[:,1],'k--',label = r'$dx_\mathrm{aft}$')
        # plt.ylim(0,1500*max(max(arc[:,0]),max(arc[:,1])))
        # plt.xlim(0,1000*(ySect[1]-ySect[0]))
        # # plt.ylim(None,1.5*max(IxSpars_vec))
        # # plt.xticks(np.array((ySect[0],yMid,ySect[1])))
        # # plt.legend(fontsize=9,ncol=3,loc=3)
        # plt.legend(fontsize=9,bbox_to_anchor=(0., 0.75, 1., .10), loc=3,
        #    ncol=4, mode="expand", borderaxespad=0.)
        # # plt.title('Spar Spacing')
        
        
        # plt.subplot(2,1,2)
        # # CS = plt.tricontourf(curvePts[:,0], curvePts[:,1]-chordwiseShift, curvePts[:,2],[0,1,2,3,5,10], cmap="YlOrRd_r")
        # # plt.colorbar(CS, orientation = "horizontal")
        # for j in range(n):
        #     for i in range(len(xSpars[j])):
        #         plt.plot(1000*yDist[j],  1000*(xSpars[j][i]-chordwiseShift),'k.',markersize=0.75)
        # # for i in [testSection-1,testSection]:
        # #     cRib = np.interp(sectBreaks[i],y_vec,np.ndarray.flatten(c_vec))
        # #     plt.plot([sectBreaks[i],sectBreaks[i]],  [-xcLL *cRib,(1-xcLL)*cRib]-chordwiseShift,'k-',linewidth=1)
        # plt.plot([-0,1000*0.5*b],  [-1000*chordwiseShift,-1000*chordwiseShift],'k-',linewidth=1)
        # plt.plot(1000*yDist,  1000*(-xcLL *c_vec-chordwiseShift),'k-',linewidth=2)
        # plt.plot(1000*yDist,1000*((1-xcLL)*c_vec-chordwiseShift),'k-',linewidth=2)
        # # plt.xticks(np.array((-0.50, -0.25, 0.0, 0.25, 0.50)))
        # plt.gca().invert_yaxis()
        # plt.gca().set_aspect(0.75)
        # plt.xlim(0,1000*(ySect[1]-ySect[0]))
        # plt.xlabel(r'$y$(mm)', fontsize=10)
        # plt.ylabel(r'$x$(mm)')
        # # plt.title('Upper Surface Curvature')
        
        if 0: # Spar Spacing
            SparSpacingFile = open("SparSpacing.txt", "w")
            j = -1
            notDone = 1
            while notDone:
                j += 1 
                ySparSpacing = 1000*yDist[j]
                if ySparSpacing >=0:
                    if ySparSpacing >1000*(ySect[1]-ySect[0]):
                        break
                    else:    
                        SparSpacingFile.write(str(ySparSpacing) +','+str(1000*arc[j,0]) +','+str(1000*arc[j,1])+'\n')
            SparSpacingFile.close()
            
        if 0: # Chord
            ChordFile = open("Chord.txt", "w")
            j = -1
            notDone = 1
            while notDone:
                j += 1 
                ySparSpacing = 1000*yDist[j]
                if ySparSpacing >=0:
                    if ySparSpacing >1000*(ySect[1]-ySect[0]):
                        break
                    else:    
                        ChordFile.write(str(ySparSpacing) + ',' + str(1000*c_vec[j,0]) + '\n')
            SparSpacingFile.close()
    if 1:
        import matplotlib.pyplot as plt
        # Generate Plot of Structure that this function will Return
        IxReq_vecPlot = np.array(IxReq_vec)*1000**4
        IxSpars_vecPlot = np.array(IxSpars_vec)*1000**4
        IxSkinPlot = np.array(IxSkin)*1000**4
        
        chordwiseShift = (-xcLL *c_vec)[n//2,0]
        width  = 5
        height = 7
        figReturn = plt.figure(figsize=(width, height))
        plt.subplot(411)
        plt.grid()
        plt.ylabel(r'Spacing Factor', fontsize=10)
        plt.xlim(0,0.5)
        plt.plot(y_vec, m_vec,'k-',label = r'$m$')
        plt.xticks(np.array((0,0.25,0.50)), ())
        plt.ylim(None,1.2*max(m_vec))
        plt.title('b = %6.3f m' %b)
        
        plt.subplot(412)
        plt.grid()
        for j in range(n):
            for i in range(len(xSparsI[j])):
                plt.plot(y_vec[j],  xSparsI[j][i]-chordwiseShift,'k.',markersize=0.75)
        for i in range(len(sectBreaks)):
            cRib = np.interp(sectBreaks[i],y_vec,np.ndarray.flatten(c_vec))
            plt.plot([sectBreaks[i],sectBreaks[i]],  [-xcLL *cRib,(1-xcLL)*cRib]-chordwiseShift,'k-',linewidth=1)
        plt.plot([-0.5,0.5],  [-chordwiseShift,-chordwiseShift],'k-',linewidth=1)
        plt.plot(y_vec,  -xcLL *c_vec-chordwiseShift,'k-',linewidth=2)
        plt.plot(y_vec,(1-xcLL)*c_vec-chordwiseShift,'k-',linewidth=2)
        plt.xlim(0,0.5)
        plt.xticks(np.array((0.0, 0.25, 0.50)))
        plt.gca().invert_yaxis()
        plt.gca().set_aspect(1/b)
        # plt.xlabel(r'Spanwise Location, [y/b]', fontsize=10)
        plt.ylabel('Chordwise \nLocation\n'+r'[x/c]')
        
        
        plt.subplot(413)
        plt.grid()
        plt.ylabel(r'$I_{xx}$ , $\mathrm{mm}^4$', fontsize=10)
        plt.xlim(0,0.5)
        plt.plot(y_vec, IxReq_vecPlot  ,'-',color='r',label = r'Req.')
        plt.plot(y_vec, IxSkinPlot     ,'-.',color=(0.9290, 0.6940, 0.1250),label = r'Skin')
        plt.plot(y_vec, IxSpars_vecPlot,'--' ,color='b',label = r'Spars')
        plt.ylim(None,1.5*max(IxSkinPlot))
        plt.xticks(np.array((0,0.25,0.50)), ())
        # plt.legend(fontsize=9,ncol=3,loc=3)
        plt.legend(fontsize=9,bbox_to_anchor=(0., 0.75, 1., .10), loc=3,
           ncol=4, mode="expand", borderaxespad=0.)
        # plt.xlabel(r'Spanwise Location, $y/b$', fontsize=10)
        
        plt.subplot(414)
        plt.grid()
        plt.ylabel(r'Factor of Safety', fontsize=10)
        plt.xlim(0,0.5)
        plt.plot(y_vec, FS_Tot      ,'k-',label = r'Total')
        plt.plot(y_vec, FS_Skin_vec ,'-.',color=(0.9290, 0.6940, 0.1250),label = r'Skin')
        plt.plot(y_vec, FS_Spars_vec,'b--',label = r'Spars')
        plt.plot([0,0.5], [FS,FS]      ,'r:',label = r'Req.')
        plt.xticks(np.array((0,0.25,0.50)))
        plt.xlabel(r'Spanwise Location, $y/b$', fontsize=10)
        plt.ylim(0,2*FS_Tot[n//2])
        # plt.legend(fontsize=9)
        plt.legend(fontsize=9,bbox_to_anchor=(0., 0.9, 1., .10), loc=2,
           ncol=3, mode="expand", borderaxespad=0.)
        plt.close(figReturn)
    
    return [xcSpars,wStruct_vec,IxTot,FSroot, Ixroot,figReturn]
    # return [xcSparsI,wStruct_vec,IxTot,FSroot, Ixroot,figReturn]

def SimpleSpars(WingShape,mb_vec,AFcoords,b,xcLL,sigma,rho_m,g,w_e):
    
    #import a whole buncha stuff
    import numpy as np
    import matplotlib.pyplot as plt
    from CodeFiles.Aerodynamics.Airfoils.operations import AFGeomAnalysis
    from CodeFiles.Structures.functions import sectionSparPlace
    from CodeFiles.Structures.functions import wingSparNormPlace
    from CodeFiles.Structures.functions import wingSparAbsPlace
    from CodeFiles.Structures.functions import rSpars
    from CodeFiles.Structures.functions import bendingInertia
    from CodeFiles.Structures.functions import sectionBendingInertia
    from CodeFiles.Structures.functions import bisection
    
    n = len(WingShape[:,0])
    y_vec     = WingShape[:,0]
    c_vec = np.zeros((n,1))
    C     = np.zeros((n,n)) 
    for j in range(n):
        c_vec[j,0]     = WingShape[j,1]
        C[j,j]     = WingShape[j,1]
    
    # Define Thickness function of AF
    AFGeom = AFGeomAnalysis(AFcoords,101)
    tauDist = np.zeros((101,2))
    tauDist[:,0] = AFGeom[:,0]
    tauDist[:,1] = AFGeom[:,3]
    taumax  = max(tauDist[:,1])
    xtmax   = tauDist[np.argmax(tauDist[:,1]),0]
        
    # IxReq_vec = Ireq_f_mb(c_vec,mb_vec,taumax,sigma)  
    IxReq_vec = (1/(2*sigma))*taumax*C@mb_vec
    
    # Required functions for optimizer
    # costFunc    = lambda m_vec: (1/n)*sum(rSpars(c_vec,wingSparNormPlace(m_vec,xtmax),tauDist))[0]
    IxSparsLam  = lambda m_vec: sum(bendingInertia(c_vec,wingSparNormPlace(m_vec,xtmax),tauDist,w_e).T)
    # IxConstr    = lambda m_vec: IxSparsLam(m_vec) - np.ndarray.flatten(IxReq_vec)
    
    import copy 
    # Initial Guess of Structure 
    IxMainSpar = (1/12)*w_e*taumax**3*c_vec**3
    if (IxMainSpar>IxReq_vec).all():
        m_vec0 = np.zeros((n,1))
        print("MAIN SPAR STRONG ENOUGH STRUCTURE")
    else:
        IxReq_max = np.max(IxReq_vec)
        iMaxStruct = np.argmax(IxReq_vec)
        maxStructSolve = lambda m: 10**10*(sectionBendingInertia(c_vec[iMaxStruct,0],sectionSparPlace(m,xtmax),tauDist,w_e) - IxReq_max)
        mMax = bisection(maxStructSolve,0,20,20)
        print('mMax = ',mMax)
        m_vec0 = copy.deepcopy(mb_vec)
        if any(IxMainSpar > IxReq_vec):
            m_vec0[IxMainSpar > IxReq_vec] = max(m_vec0[IxMainSpar > IxReq_vec])
        m_vec0 = m_vec0 - min(m_vec0)
        m_vec0 = mMax*m_vec0/max(m_vec0)
    
    m_vec       = copy.deepcopy(m_vec0)
    xcSpars     = wingSparNormPlace(m_vec,xtmax)
    xSpars   = wingSparAbsPlace(xcSpars,c_vec,xcLL)
    IxSpars_vec = IxSparsLam(m_vec)
    
    # Final Tweeks on Structure
    for j in range(n):  #Sweep from left tip to right tip
        if IxReq_vec[j,0]>IxSpars_vec[j]:
            maxStructSolve = lambda m: 10**10*(sectionBendingInertia(c_vec[j,0],sectionSparPlace(m,xtmax),tauDist,w_e) - IxReq_vec[j,0])
            m_vec[j,0] = bisection(maxStructSolve,0,20,20)
    
    xcSpars     = wingSparNormPlace(m_vec,xtmax)
    xSpars   = wingSparAbsPlace(xcSpars,c_vec,xcLL)
    IxSpars_vec = IxSparsLam(m_vec)
    rSpars_vec = rSpars(c_vec,xcSpars,tauDist)
    wStruct_vec = rho_m*g*w_e*rSpars_vec
    
    width  = 3.25
    height = 3.25
    plt.figure(figsize=(width, height))
    
    plt.subplot(211)
    plt.grid()
    plt.ylabel(r'$I_x$', fontsize=10)
    plt.xlim(-0.5,0.5)
    plt.plot(y_vec, IxReq_vec  ,'k--',label = r'$I_\mathrm{req}$')
    plt.plot(y_vec, IxSpars_vec,'k-',label = r'$I_\mathrm{Struct}$')
    plt.plot(y_vec, IxMainSpar ,'k:',label = r'$I_\mathrm{MainSpar}$')
    plt.ylim(None,1.5*max(IxSpars_vec))
    plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)), ())
    # plt.legend(fontsize=9,ncol=3,loc=3)
    plt.legend(fontsize=9,bbox_to_anchor=(0., 0.75, 1., .10), loc=3,
       ncol=3, mode="expand", borderaxespad=0.)
    
    plt.subplot(212)
    plt.grid()
    plt.ylabel(r'm', fontsize=10)
    plt.xlim(-0.5,0.5)
    plt.plot(y_vec, m_vec0,'k--',label = r'$m_\mathrm{0}$')
    plt.plot(y_vec, m_vec,'k-',label = r'$m$')
    plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)))
    plt.xlabel(r'Spanwise Location, $y/b$', fontsize=10)
    plt.ylim(None,1.5*max(m_vec))
    # plt.legend(fontsize=9)
    plt.legend(fontsize=9,bbox_to_anchor=(0., 0.9, 1., .10), loc=2,
       ncol=3, mode="expand", borderaxespad=0.)
         
    width  = 3.25
    height = 3
    plt.figure(figsize=(width, height))
    
    plt.grid()
    plt.ylabel(r'$z/c$', fontsize=10)
    plt.xlabel(r'$x/c$', fontsize=10)
    plt.axis('equal')
    plt.plot(AFGeom[:,0], AFGeom[:,1], linestyle='-')
    plt.plot(AFGeom[:,0], AFGeom[:,2], linestyle='-')
    plt.plot(AFGeom[:,0], AFGeom[:,4], linestyle='--')
    plt.xticks(np.array((0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)))
    # plt.yticks(np.array((0,0.05,0.10)))
    # plt.ylim(0, 1.2*max(c_vec))
    plt.xlim(0.0,1.0)
    
    chordwiseShift = (-xcLL *c_vec)[n//2,0]
    width  = 3.25
    height = 3.5
    plt.figure(figsize=(width, height))
    plt.subplot(211)
    plt.grid()
    for i in range(len(xcSpars[0,:])):
        plt.plot(y_vec,  xSpars[:,i]-chordwiseShift)
    plt.plot(y_vec,  -xcLL *c_vec-chordwiseShift,'k-')
    plt.plot(y_vec,(1-xcLL)*c_vec-chordwiseShift,'k-')
    plt.xlim(0,0.5)
    plt.gca().invert_yaxis()
    plt.gca().set_aspect("equal")
    plt.xlim(0.0,0.5)
    plt.xlabel(r'Spanwise Location, $y/b$', fontsize=10)
    plt.ylabel(r'$c/b$')
    
    plt.subplot(212)
    plt.grid()
    for i in range(len(xcSpars[0,:])):
        plt.plot(y_vec,  xSpars[:,i]-chordwiseShift)
    plt.plot(y_vec,  (-xcLL *c_vec)-chordwiseShift,'k-')
    plt.plot(y_vec,((1-xcLL)*c_vec)-chordwiseShift,'k-')
    plt.xlim(-0.5,0.0)
    plt.gca().invert_yaxis()
    plt.gca().set_aspect("equal")
    plt.xlim(-0.5,0.0)
    plt.xlabel(r'Spanwise Location, $y/b$', fontsize=10)
    plt.ylabel(r'$c/b$')
    
    return [xcSpars,wStruct_vec]

def SimpleSparsWAilSkinInDesign(WingShape,mb_vec,AFcoords,b,xcLL,sigma,rho_m,g,w_e):
    
    #import a whole buncha stuff
    import numpy as np
    import matplotlib.pyplot as plt
    from CodeFiles.Aerodynamics.Airfoils.operations import AFGeomAnalysis
    from CodeFiles.Aerodynamics.Airfoils.operations import AFBendProperties
    from CodeFiles.Structures.functions import sectionSparPlace
    from CodeFiles.Structures.functions import wingSparNormPlace
    from CodeFiles.Structures.functions import wingSparAbsPlace
    from CodeFiles.Structures.functions import rSpars
    from CodeFiles.Structures.functions import bendingInertia
    from CodeFiles.Structures.functions import sectionBendingInertia
    from CodeFiles.Structures.functions import bisection
    
    SkinStrength = 0.8
    
    [IxAFnorm,zcbar] = AFBendProperties(AFcoords,xcLL)
    n = len(WingShape[:,0])
    y_vec     = WingShape[:,0]
    c_vec = np.zeros((n,1))
    C     = np.zeros((n,n))
    IxSkin = np.zeros((n,1))
    for j in range(n):
        c_vec[j,0]  = b*WingShape[j,1] # NOTE c is now DIMENSIONAL
        C[j,j]      = b*WingShape[j,1]
        IxSkin[j,0] = SkinStrength*w_e*(c_vec[j,0])**3*IxAFnorm #skin bending intertia
    
    # Define Thickness function of AF
    AFGeom = AFGeomAnalysis(AFcoords,101)
    tauDist = np.zeros((101,2))
    tauDist[:,0] = AFGeom[:,0]
    tauDist[:,1] = AFGeom[:,3]
    taumax  = max(tauDist[:,1])
    xtmax   = tauDist[np.argmax(tauDist[:,1]),0]
        
    # IxReq_vec = Ireq_f_mb(c_vec,mb_vec,taumax,sigma)  
    IxReq_vec = (1/(2*sigma))*taumax*C@mb_vec
    IxReqSpars_vec = IxReq_vec - IxSkin
    
    # Required functions for optimizer
    # costFunc    = lambda m_vec: (1/n)*sum(rSpars(c_vec,wingSparNormPlace(m_vec,xtmax),tauDist))[0]
    IxSparsLam  = lambda m_vec: sum(bendingInertia(c_vec,wingSparNormPlace(m_vec,xtmax),tauDist,w_e).T)
    # IxConstr    = lambda m_vec: IxSparsLam(m_vec) - np.ndarray.flatten(IxReq_vec)
    
    import copy 
    # Initial Guess of Structure
    IxMainSpar = (1/12)*w_e*(taumax*c_vec)**3
    if (IxReqSpars_vec>IxReq_vec).all():
        m_vec0 = np.zeros((n,1))
        print("SKIN + MAIN SPAR STRONG ENOUGH STRUCTURE")
    else:
        IxSparsReq_max = np.max(IxReqSpars_vec)
        iMaxStruct = np.argmax(IxReqSpars_vec)
        # breakpoint()
        maxStructSolve = lambda m: 10**10*(sectionBendingInertia(c_vec[iMaxStruct,0],sectionSparPlace(m,xtmax),tauDist,w_e) - IxSparsReq_max)
        mMax = bisection(maxStructSolve,0,40,20)
        if mMax == None:
            print("STRUCTURE NOT POSSIBLE, LOAD TOO HIGH!!!!!!!")
        elif mMax == 0:
            m_vec0 = np.zeros((n,1))
        else:
            print('mMax = ',mMax)
            m_vec0 = copy.deepcopy(mb_vec)
            if any(IxMainSpar > IxReqSpars_vec):
                m_vec0[IxMainSpar > IxReqSpars_vec] = max(m_vec0[IxMainSpar > IxReq_vec])
            m_vec0 = m_vec0 - min(m_vec0)
            # breakpoint()
            m_vec0 = mMax*m_vec0/max(m_vec0)
    
    m_vec       = copy.deepcopy(m_vec0)
    xcSpars     = wingSparNormPlace(m_vec,xtmax)
    xSpars   = wingSparAbsPlace(xcSpars,c_vec,xcLL)
    IxSpars_vec = IxSparsLam(m_vec)
    
    # Final Tweeks on Structure
    for j in range(n):  #Sweep from left tip to right tip
        if IxReqSpars_vec[j,0]>IxSpars_vec[j]:
            maxStructSolve = lambda m: 10**10*(sectionBendingInertia(c_vec[j,0],sectionSparPlace(m,xtmax),tauDist,w_e) - IxReqSpars_vec[j,0])
            m_vec[j,0] = bisection(maxStructSolve,0,20,20)
    
    xcSpars     = wingSparNormPlace(m_vec,xtmax)
    xSpars   = wingSparAbsPlace(xcSpars,c_vec,xcLL)
    IxSpars_vec = IxSparsLam(m_vec)
    rSpars_vec = rSpars(c_vec,xcSpars,tauDist)
    wStruct_vec = rho_m*g*w_e*rSpars_vec
    
    
    width  = 3.25
    height = 3.25
    plt.figure(figsize=(width, height))
    
    plt.subplot(211)
    plt.grid()
    plt.ylabel(r'$I_x$', fontsize=10)
    plt.xlim(-0.5,0.5)
    plt.plot(y_vec, IxReq_vec  ,'k--',label = r'$I_\mathrm{req}$')
    plt.plot(y_vec, IxSpars_vec,'k-',label = r'$I_\mathrm{Spars}$')
    plt.plot(y_vec, IxMainSpar ,'k:',label = r'$I_\mathrm{MainSpar}$')
    plt.plot(y_vec, IxSkin     ,'r-.',label = r'$I_\mathrm{Skin}$')
    # plt.ylim(None,1.5*max(IxSpars_vec))
    plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)), ())
    # plt.legend(fontsize=9,ncol=3,loc=3)
    plt.legend(fontsize=9,bbox_to_anchor=(0., 0.75, 1., .10), loc=3,
       ncol=4, mode="expand", borderaxespad=0.)
    
    plt.subplot(212)
    plt.grid()
    plt.ylabel(r'm', fontsize=10)
    plt.xlim(-0.5,0.5)
    plt.plot(y_vec, m_vec0,'k--',label = r'$m_\mathrm{0}$')
    plt.plot(y_vec, m_vec,'k-',label = r'$m$')
    plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)))
    plt.xlabel(r'Spanwise Location, $y/b$', fontsize=10)
    plt.ylim(None,1.5*max(m_vec))
    # plt.legend(fontsize=9)
    plt.legend(fontsize=9,bbox_to_anchor=(0., 0.9, 1., .10), loc=2,
       ncol=3, mode="expand", borderaxespad=0.)
         
    width  = 3.25
    height = 3
    plt.figure(figsize=(width, height))
    
    plt.grid()
    plt.ylabel(r'$z/c$', fontsize=10)
    plt.xlabel(r'$x/c$', fontsize=10)
    plt.axis('equal')
    plt.plot(AFGeom[:,0], AFGeom[:,1], linestyle='-')
    plt.plot(AFGeom[:,0], AFGeom[:,2], linestyle='-')
    plt.plot(AFGeom[:,0], AFGeom[:,4], linestyle='--')
    plt.xticks(np.array((0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)))
    # plt.yticks(np.array((0,0.05,0.10)))
    # plt.ylim(0, 1.2*max(c_vec))
    plt.xlim(0.0,1.0)
    
    chordwiseShift = (-xcLL *c_vec)[n//2,0]
    width  = 3.25
    height = 3.5
    plt.figure(figsize=(width, height))
    plt.subplot(211)
    plt.grid()
    for i in range(len(xcSpars[0,:])):
        plt.plot(y_vec,  xSpars[:,i]-chordwiseShift)
    plt.plot(y_vec,  -xcLL *c_vec-chordwiseShift,'k-')
    plt.plot(y_vec,(1-xcLL)*c_vec-chordwiseShift,'k-')
    plt.xlim(0,0.5)
    plt.gca().invert_yaxis()
    plt.gca().set_aspect("equal")
    plt.xlim(0.0,0.5)
    plt.xlabel(r'Spanwise Location, $y/b$', fontsize=10)
    plt.ylabel(r'$c/b$')
    
    plt.subplot(212)
    plt.grid()
    for i in range(len(xcSpars[0,:])):
        plt.plot(y_vec,  xSpars[:,i]-chordwiseShift)
    plt.plot(y_vec,  (-xcLL *c_vec)-chordwiseShift,'k-')
    plt.plot(y_vec,((1-xcLL)*c_vec)-chordwiseShift,'k-')
    plt.xlim(-0.5,0.0)
    plt.gca().invert_yaxis()
    plt.gca().set_aspect("equal")
    plt.xlim(-0.5,0.0)
    plt.xlabel(r'Spanwise Location, $y/b$', fontsize=10)
    plt.ylabel(r'$c/b$')
    
    return [xcSpars,wStruct_vec]

def SimpleSparsWAilNoSkin(WingShape,mb_vec,AFcoords,b,xcLL,sigma,rho_m,g,w_e):
    
        
    
    #import a whole buncha stuff
    import numpy as np
    import matplotlib.pyplot as plt
    from CodeFiles.Aerodynamics.Airfoils.operations import AFGeomAnalysis
    from CodeFiles.Structures.functions import sectionSparPlace
    from CodeFiles.Structures.functions import wingSparNormPlace
    from CodeFiles.Structures.functions import wingSparAbsPlace
    from CodeFiles.Structures.functions import rSpars
    from CodeFiles.Structures.functions import bendingInertia
    from CodeFiles.Structures.functions import sectionBendingInertia
    from CodeFiles.Structures.functions import bisection
    
    n = len(WingShape[:,0])
    y_vec     = WingShape[:,0]
    c_vec = np.zeros((n,1))
    C     = np.zeros((n,n)) 
    for j in range(n):
        c_vec[j,0]     = WingShape[j,1]
        C[j,j]     = WingShape[j,1]
    
    # Define Thickness function of AF
    AFGeom = AFGeomAnalysis(AFcoords,101)
    tauDist = np.zeros((101,2))
    tauDist[:,0] = AFGeom[:,0]
    tauDist[:,1] = AFGeom[:,3]
    taumax  = max(tauDist[:,1])
    xtmax   = tauDist[np.argmax(tauDist[:,1]),0]
        
    # IxReq_vec = Ireq_f_mb(c_vec,mb_vec,taumax,sigma)  
    IxReq_vec = (1/(2*sigma))*taumax*C@mb_vec
    
    # Required functions for optimizer
    # costFunc    = lambda m_vec: (1/n)*sum(rSpars(c_vec,wingSparNormPlace(m_vec,xtmax),tauDist))[0]
    IxSparsLam  = lambda m_vec: sum(bendingInertia(c_vec,wingSparNormPlace(m_vec,xtmax),tauDist,w_e).T)
    # IxConstr    = lambda m_vec: IxSparsLam(m_vec) - np.ndarray.flatten(IxReq_vec)
    
    import copy 
    # Initial Guess of Structure 
    IxMainSpar = (1/12)*w_e*taumax**3*c_vec**3
    if (IxMainSpar>IxReq_vec).all():
        m_vec0 = np.zeros((n,1))
        print("MAIN SPAR STRONG ENOUGH STRUCTURE")
    else:
        IxReq_max = np.max(IxReq_vec)
        iMaxStruct = np.argmax(IxReq_vec)
        maxStructSolve = lambda m: 10**10*(sectionBendingInertia(c_vec[iMaxStruct,0],sectionSparPlace(m,xtmax),tauDist,w_e) - IxReq_max)
        mMax = bisection(maxStructSolve,0,20,20)
        print('mMax = ',mMax)
        m_vec0 = copy.deepcopy(mb_vec)
        if any(IxMainSpar > IxReq_vec):
            m_vec0[IxMainSpar > IxReq_vec] = max(m_vec0[IxMainSpar > IxReq_vec])
        m_vec0 = m_vec0 - min(m_vec0)
        m_vec0 = mMax*m_vec0/max(m_vec0)
    
    m_vec       = copy.deepcopy(m_vec0)
    xcSpars     = wingSparNormPlace(m_vec,xtmax)
    xSpars   = wingSparAbsPlace(xcSpars,c_vec,xcLL)
    IxSpars_vec = IxSparsLam(m_vec)
    
    # Final Tweeks on Structure
    for j in range(n):  #Sweep from left tip to right tip
        if IxReq_vec[j,0]>IxSpars_vec[j]:
            maxStructSolve = lambda m: 10**10*(sectionBendingInertia(c_vec[j,0],sectionSparPlace(m,xtmax),tauDist,w_e) - IxReq_vec[j,0])
            m_vec[j,0] = bisection(maxStructSolve,0,20,20)
    
    xcSpars     = wingSparNormPlace(m_vec,xtmax)
    xSpars   = wingSparAbsPlace(xcSpars,c_vec,xcLL)
    IxSpars_vec = IxSparsLam(m_vec)
    rSpars_vec = rSpars(c_vec,xcSpars,tauDist)
    wStruct_vec = rho_m*g*w_e*rSpars_vec
    
    
    width  = 3.25
    height = 3.25
    plt.figure(figsize=(width, height))
    
    plt.subplot(211)
    plt.grid()
    plt.ylabel(r'$I_x$', fontsize=10)
    plt.xlim(-0.5,0.5)
    plt.plot(y_vec, IxReq_vec  ,'k--',label = r'$I_\mathrm{req}$')
    plt.plot(y_vec, IxSpars_vec,'k-',label = r'$I_\mathrm{Struct}$')
    plt.plot(y_vec, IxMainSpar ,'k:',label = r'$I_\mathrm{MainSpar}$')
    plt.ylim(None,1.5*max(IxSpars_vec))
    plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)), ())
    # plt.legend(fontsize=9,ncol=3,loc=3)
    plt.legend(fontsize=9,bbox_to_anchor=(0., 0.75, 1., .10), loc=3,
       ncol=3, mode="expand", borderaxespad=0.)
    
    plt.subplot(212)
    plt.grid()
    plt.ylabel(r'm', fontsize=10)
    plt.xlim(-0.5,0.5)
    plt.plot(y_vec, m_vec0,'k--',label = r'$m_\mathrm{0}$')
    plt.plot(y_vec, m_vec,'k-',label = r'$m$')
    plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)))
    plt.xlabel(r'Spanwise Location, $y/b$', fontsize=10)
    plt.ylim(None,1.5*max(m_vec))
    # plt.legend(fontsize=9)
    plt.legend(fontsize=9,bbox_to_anchor=(0., 0.9, 1., .10), loc=2,
       ncol=3, mode="expand", borderaxespad=0.)
         
    width  = 3.25
    height = 3
    plt.figure(figsize=(width, height))
    
    plt.grid()
    plt.ylabel(r'$z/c$', fontsize=10)
    plt.xlabel(r'$x/c$', fontsize=10)
    plt.axis('equal')
    plt.plot(AFGeom[:,0], AFGeom[:,1], linestyle='-')
    plt.plot(AFGeom[:,0], AFGeom[:,2], linestyle='-')
    plt.plot(AFGeom[:,0], AFGeom[:,4], linestyle='--')
    plt.xticks(np.array((0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)))
    # plt.yticks(np.array((0,0.05,0.10)))
    # plt.ylim(0, 1.2*max(c_vec))
    plt.xlim(0.0,1.0)
    
    chordwiseShift = (-xcLL *c_vec)[n//2,0]
    width  = 3.25
    height = 3.5
    plt.figure(figsize=(width, height))
    plt.subplot(211)
    plt.grid()
    for i in range(len(xcSpars[0,:])):
        plt.plot(y_vec,  xSpars[:,i]-chordwiseShift)
    plt.plot(y_vec,  -xcLL *c_vec-chordwiseShift,'k-')
    plt.plot(y_vec,(1-xcLL)*c_vec-chordwiseShift,'k-')
    plt.xlim(0,0.5)
    plt.gca().invert_yaxis()
    plt.gca().set_aspect("equal")
    plt.xlim(0.0,0.5)
    plt.xlabel(r'Spanwise Location, $y/b$', fontsize=10)
    plt.ylabel(r'$c/b$')
    
    plt.subplot(212)
    plt.grid()
    for i in range(len(xcSpars[0,:])):
        plt.plot(y_vec,  xSpars[:,i]-chordwiseShift)
    plt.plot(y_vec,  (-xcLL *c_vec)-chordwiseShift,'k-')
    plt.plot(y_vec,((1-xcLL)*c_vec)-chordwiseShift,'k-')
    plt.xlim(-0.5,0.0)
    plt.gca().invert_yaxis()
    plt.gca().set_aspect("equal")
    plt.xlim(-0.5,0.0)
    plt.xlabel(r'Spanwise Location, $y/b$', fontsize=10)
    plt.ylabel(r'$c/b$')
    
    return [xcSpars,wStruct_vec]