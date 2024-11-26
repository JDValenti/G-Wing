def PlotWing(y_vec,c_vec,alpha_vec,b,xcLL):
    
    import numpy as np
    import matplotlib.pyplot as plt
    
    n = len(y_vec)
    chordwiseShift = (-xcLL *c_vec)[n//2]
    print("xcLL        = ",xcLL)
    print("c_vec[n//2] = ",c_vec[n//2])
    print("Chordshift  = ",chordwiseShift)
    
    width  = 6.5
    height = 6
    plt.figure(figsize=(width, height))
    
    plt.subplot(311)
    plt.grid()
    plt.ylabel(r'Chord, $c/b$', fontsize=18)
    # plt.axis('equal')
    # plt.plot(y, cvec[:,0],'k:')
    plt.plot(y_vec, b*c_vec,'k-')
    plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)), ())
    plt.yticks(fontsize=18)
    # plt.yticks(np.array((0,0.05,0.10)))
    # plt.ylim(0, 1.2*max(cvec[:,-1]))
    plt.xlim(-0.5,0.5)
    
    
    plt.subplot(312)
    plt.grid()
    plt.ylabel(r'Twist (deg)', fontsize=18)
    # plt.ylim(0, 1.2*(180/pi)*max(alpha))
    plt.xlim(-0.5,0.5)
    plt.plot(y_vec, (180/np.pi)*(alpha_vec - alpha_vec[n//2]),'k-')
    plt.xticks(np.array((-0.50,-0.25,0,0.25,0.50)), ())
    plt.yticks(fontsize=18)
    

    plt.subplot(313)
    # plt.grid()
    plt.plot([-0.5,0.5],  [b*-chordwiseShift,b*-chordwiseShift],'k--',linewidth=2)
    plt.plot(y_vec,b*(  -xcLL *c_vec-chordwiseShift),'k-',linewidth=2)
    plt.plot(y_vec,b*((1-xcLL)*c_vec-chordwiseShift),'k-',linewidth=2)
    plt.xlim(0,0.5)
    plt.xticks(np.array((-0.50, -0.25, 0.0, 0.25, 0.50)))
    plt.gca().invert_yaxis()
    plt.gca().set_aspect(1/b)
    plt.xticks(np.array((-0.5, -0.25, 0.0, 0.25, 0.5)),fontsize=18)
    plt.yticks(np.array((0.0, 0.075, 0.15)),fontsize=18)
    plt.xlabel(r'Spanwise Location, $y/b$', fontsize=20)
    plt.ylabel(r'$x/b$', fontsize=20)
        
def LLanalyzeAlpha(WingShape,A,alpha_w,dispPlots):
    """
    Inports:
        * WingShape
            type: n x 3 numpy array
            description: defines wing shape
                rows correspond to spanwise location
                columns correspond to engineering values:
                    col  0: y
                    col  1: c/b
                    col  2: alpha [rad]
        * A
            type:  float or matrix
            description: defines sectional lift curve slope.  If float, then 
                constant value for entire span. If matrix, then spanwise 
                distribution
        * alpha_w
            type: float or array
            description: alpha_w values at which to evaluate wing
        * disp plot
            type: integer (0 or 1)
            description: turns plotting on (1) or off (0)
            
    Returns:
        * vectors
            type:  n x 3 x ncL numpy array
            description: contain the spanwise distributions of 
                (col 0) nonDim Vorticity, (col 1) cl, and (col 2) alpha_i for 
                each cL condition
        * scalars
            type:  3 x ncL numpy array
            description: contain the scalar values of (row 0) cL, (row 1) cDi, 
                and (row 2) span efficiency for each cL condition
    """
    import numpy as np
    from numpy.linalg import inv
    from numpy import pi
    
    print("Lifting Line Analysis----------------------")
    n = len(WingShape)
    # Define needed vectors and matrices ------------------------
    
    WingShape   = np.array(WingShape)
    y_temp       = WingShape[:,0]
    c_temp       = WingShape[:,1]
    alpha_r     = WingShape[n//2,2]
    alphatw_temp = WingShape[:,2] - alpha_r
    
    y_vec = np.zeros((n,1))
    c_vec = np.zeros((n,1))
    alphatw_vec  = np.zeros((n,1))
    for i in range(n):
        y_vec[i,0] = y_temp[i]
        c_vec[i,0] = c_temp[i]
        alphatw_vec[i,0] = alphatw_temp[i]
    del(y_temp, c_temp, alphatw_temp)    
    
    C       = np.diag(np.ndarray.flatten(c_vec))
    Cinv    = np.diag(np.ndarray.flatten(1/c_vec))
    
    e = np.ones((n,1))
    I = np.eye(n)
    
    Q = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            Q[i,j] = 1/(1-4*(i-j)**2)
    Qinv = inv(Q)
    print("Q-Matrix Generated!\n")
    if np.shape(A) == ():
        W = ((A*n)/(2*pi))*Q
        A = A*I
    else:
        W = ((n)/(2*pi))*A@Q
        
    CinvWinv = inv(Cinv + W)
    eCe = e.T@C@e
    AR = n/eCe
    
    alphaSeq  = np.array(alpha_w)
    if np.size(alphaSeq) == 1:
        nSeq = 1
    else:
        nSeq      = len(alphaSeq)
    cL        = 0
    cDi       = 0
    e_span    = 0
    Gamma_vec = np.zeros((n,1))
    cl_vec    = np.zeros((n,1))
    w_vec     = np.zeros((n,1))
    
    vectors   = np.zeros((n,3,nSeq))
    scalars   = np.zeros((3,nSeq))
    
    for i in range(nSeq):
        if nSeq == 1:
            Gamma_vec = 0.5*(CinvWinv)@A@(alphatw_vec + alphaSeq*e)
        else:
            Gamma_vec = 0.5*(CinvWinv)@A@(alphatw_vec + alphaSeq[i]*e)
        cl_vec    = 2*Cinv@Gamma_vec
        w_vec     = (n/np.pi)*Q@Gamma_vec
        
        vectors[:,0,i] = Gamma_vec[:,0]
        vectors[:,1,i] = cl_vec[:,0]
        vectors[:,2,i] = w_vec[:,0]
        
        cL  = (2*e.T@Gamma_vec/eCe)[0,0]
        cDi = (((2*n)/(np.pi*eCe))*Gamma_vec.T@Q@Gamma_vec)[0,0]
        e_span = (cL**2/(np.pi*AR*cDi))[0,0]
        
        print("alpha_w = ", alphaSeq[i]*180/pi, "deg;  cL = ", cL, ";  cDi = ", cDi, ";  e = ", e_span)
        scalars[0,i] = cL
        scalars[1,i] = cDi
        scalars[2,i] = e_span
        
            
    # print("cLmax = ",cLmax)
    # l_vec_clmax = l_vec[:,kcLmax]
    
    if dispPlots:
    
        import matplotlib.pyplot as plt
        lineStyleOrder = ['--',':','-.']
        
        width  = 3.25
        height = 3.5
        plt.figure(figsize=(width, height))
        plt.subplot(311)
        plt.grid()
        plt.ylabel(r'$c c_\ell/c_\mathrm{ref}$', fontsize=10)
        # plt.axis('equal')
        for i in range(nSeq):
            plt.plot(y_vec, vectors[:,1,i]*c_vec[:,0]*AR[0,0],lineStyleOrder[i%3])
        plt.xticks(np.array((-0.5,-0.25,0,0.25,0.5)), ())
        plt.xlim([-0.5,0.5])
        plt.xlabel(r'$y/b$', fontsize=10)
        
        plt.subplot(312)
        plt.grid()
        plt.ylabel(r'$c_\ell$', fontsize=10)
        # plt.axis('equal')
        for i in range(nSeq):
            plt.plot(y_vec, vectors[:,1,i],lineStyleOrder[i%3])
        plt.xticks(np.array((-0.5,-0.25,0,0.25,0.5)), ())
        plt.xlim([-0.5,0.5])
        plt.xlabel(r'$y/b$', fontsize=10)
        
        plt.subplot(313)
        plt.grid()
        plt.ylabel(r'$\alpha_i$', fontsize=10)
        # plt.axis('equal')
        for i in range(nSeq):
            plt.plot(y_vec, vectors[:,2,i]*180/np.pi,lineStyleOrder[i%3])
        plt.xticks(np.array((-0.5,-0.25,0,0.25,0.5)))
        plt.xlim([-0.5,0.5])
        plt.xlabel(r'$y/b$', fontsize=10)
        
    return(vectors,scalars)

def LLanalyzeCL(WingShape,A,cL,dispPlots):
    """
    Inports:
        * WingShape
            type: n x 3 numpy array
            description: defines wing shape
                rows correspond to spanwise location
                columns correspond to engineering values:
                    col  0: y
                    col  1: c/b
                    col  2: alpha [rad]
        * A
            type:  float or matrix
            description: defines sectional lift curve slope.  If float, then 
                constant value for entire span. If matrix, then spanwise 
                distribution
        * cL
            type: float or array
            description: cL values at which to evaluate wing
        * disp plot
            type: integer (0 or 1)
            description: turns plotting on (1) or off (0)
            
    Returns:
        * vectors
            type:  n x 3 x ncL numpy array
            description: contain the spanwise distributions of 
                (col 0) nonDim Vorticity, (col 1) cl, and (col 2) alpha_i for 
                each cL condition
        * scalars
            type:  3 x ncL numpy array
            description: contain the scalar values of (row 0) cL, (row 1) cDi, 
                and (row 2) span efficiency for each cL condition
    """
    import numpy as np
    from numpy.linalg import inv
    from numpy import pi
    
    print("Lifting Line Analysis----------------------")
    n = len(WingShape)
    # Define needed vectors and matrices ------------------------
    
    WingShape    = np.array(WingShape)
    y_temp       = WingShape[:,0]
    c_temp       = WingShape[:,1]
    alpha_r      = WingShape[n//2,2]
    alphatw_temp = WingShape[:,2] - alpha_r
    
    y_vec = np.zeros((n,1))
    c_vec = np.zeros((n,1))
    alphatw_vec  = np.zeros((n,1))
    for i in range(n):
        y_vec[i,0] = y_temp[i]
        c_vec[i,0] = c_temp[i]
        alphatw_vec[i,0] = alphatw_temp[i]
    del(y_temp, c_temp, alphatw_temp)
    
    C       = np.diag(np.ndarray.flatten(c_vec))
    Cinv    = np.diag(np.ndarray.flatten(1/c_vec))
    
    e = np.ones((n,1))
    I = np.eye(n)
    
    Q = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            Q[i,j] = 1/(1-4*(i-j)**2)
    Qinv = inv(Q)
    print("Q-Matrix Generated!\n")
    if np.shape(A) == ():
        W = ((A*n)/(2*pi))*Q
        A = A*I
    else:
        W = ((n)/(2*pi))*A@Q
        
    CinvWinv = inv(Cinv + W)
    eCe = e.T@C@e
    AR = n/eCe
    
    # Solve the alpha_w for each cL
    Gamma0_vec = 0.5*(CinvWinv)@A@(alphatw_vec)
    cL0        = (2*e.T@Gamma0_vec/eCe)[0,0]
    cLalpha    = (e.T@CinvWinv@A@e/eCe)[0,0]
    
    if np.size(cL) == 1:
        nSeq = 1
    else:
        nSeq      = len(cL)
        
    alphaSeq = np.zeros(np.shape(cL))
    for i in range(nSeq):
        alphaSeq[i] = (cL[i] - cL0)/cLalpha    
    
    cL        = 0
    cDi       = 0
    e_span    = 0
    Gamma_vec = np.zeros((n,1))
    cl_vec    = np.zeros((n,1))
    w_vec     = np.zeros((n,1))
    
    vectors   = np.zeros((n,3,nSeq))
    scalars   = np.zeros((3,nSeq))
    
    for i in range(nSeq):
        if nSeq == 1:
            Gamma_vec = 0.5*(CinvWinv)@A@(alphatw_vec + alphaSeq*e)
        else:
            Gamma_vec = 0.5*(CinvWinv)@A@(alphatw_vec + alphaSeq[i]*e)
        cl_vec    = 2*Cinv@Gamma_vec
        w_vec     = (n/np.pi)*Q@Gamma_vec
        # w_vec     = (2*n/np.pi)*Q@Gamma_vec   # NEED TO CHECK!!!!!!!!
        
        vectors[:,0,i] = Gamma_vec[:,0]
        vectors[:,1,i] = cl_vec[:,0]
        vectors[:,2,i] = w_vec[:,0]
        
        cL  = (2*e.T@Gamma_vec/eCe)[0,0]
        cDi = (((2*n)/(np.pi*eCe))*Gamma_vec.T@Q@Gamma_vec)[0,0]
        e_span = (cL**2/(np.pi*AR*cDi))[0,0]
        
        print("cL = %3.2f; alpha_w = %4.2fdeg; cDi = %5.4f; e = %3.2f"%(cL, alphaSeq[i]*180/pi, cDi, e_span))
        # print("cL = ", cL,";  alpha_w = ", alphaSeq[i]*180/pi, "deg;  cDi = ", cDi, ";  e = ", e_span)
        scalars[0,i] = cL
        scalars[1,i] = cDi
        scalars[2,i] = e_span
        
            
    # print("cLmax = ",cLmax)
    # l_vec_clmax = l_vec[:,kcLmax]
    
    if dispPlots:
    
        import matplotlib.pyplot as plt
        lineStyleOrder = ['--',':','-.']
        
        width  = 3.25
        height = 5
        plt.figure(figsize=(width, height))
        plt.subplot(211)
        plt.grid()
        plt.ylabel(r'$c c_\ell/c_\mathrm{ref}$', fontsize=10)
        # plt.axis('equal')
        maxcCL = 0
        for i in range(nSeq):
            cCLcref = vectors[:,1,i]*c_vec[:,0]*AR[0,0]
            plt.plot(y_vec, cCLcref, lineStyleOrder[i%3],
                     label = "cL = %3.2f; cDi = %5.4f; e = %3.2f"%(scalars[0,i],scalars[1,i],scalars[2,i]))
            maxcCL = max(maxcCL,max(cCLcref))
        plt.xticks(np.array((-0.5,-0.25,0,0.25,0.5)), ())
        plt.xlim([-0.5,0.5])
        plt.ylim(top=1.8*maxcCL)
        plt.legend(fontsize=8.5, loc='best',ncol=1);
        
        plt.subplot(413)
        plt.grid()
        plt.ylabel(r'$c_\ell$', fontsize=10)
        # plt.axis('equal')
        for i in range(nSeq):
            plt.plot(y_vec, vectors[:,1,i],lineStyleOrder[i%3])
        plt.xticks(np.array((-0.5,-0.25,0,0.25,0.5)), ())
        plt.xlim([-0.5,0.5])
        plt.xlabel(r'$y/b$', fontsize=10)
        
        plt.subplot(414)
        plt.grid()
        plt.ylabel(r'$\alpha_i$', fontsize=10)
        # plt.axis('equal')
        for i in range(nSeq):
            plt.plot(y_vec, vectors[:,2,i]*180/np.pi,lineStyleOrder[i%3])
        plt.ylim([0,7])
        plt.xticks(np.array((-0.5,-0.25,0,0.25,0.5)))
        plt.xlim([-0.5,0.5])
        plt.xlabel(r'$y/b$', fontsize=10)
        
    return(vectors,scalars)

def AVLwrite(fileName,b,AR,nAVL,xcBV,WingShape,Folder = 0):
    """
    Writes WingShape to AVL file
        
    @author: justi
    """    
    import os
    import shutil
    import numpy as np
    from numpy import pi
    # from CodeFiles.OpenSCAD.GWingSCADScript import switchSettings
    # from CodeFiles.OpenSCAD.GWingSCADScript import baseScript
    
    # Some basic calculations up front
    n = len(WingShape)  
    Bref = b
    # Bref = 1
    Cref = Bref/AR
    Sref = Bref*Cref
    
    # Reduce discretization to AVL Limit
    dTheta = pi/(nAVL[0] - 1)
    yLLtip = WingShape[0,0]
    
    yInterp = np.zeros((nAVL[0],1))
    WingShapeAVL = np.zeros((nAVL[0],3))
    for i in range(nAVL[0]):
        WingShapeAVL[i,0] = -yLLtip*np.cos(pi - i*dTheta)

    
    WingShapeAVL[:,1] = np.interp(WingShapeAVL[:,0],WingShape[:,0],WingShape[:,1], left=WingShape[0,1], right=WingShape[-1,1])
    WingShapeAVL[:,2] = np.interp(WingShapeAVL[:,0],WingShape[:,0],WingShape[:,2], left=WingShape[0,2], right=WingShape[-1,2])
    
    # breakpoint()
    
    # import matplotlib.pyplot as plt
    # plt.figure()
    # plt.plot(WingShapeAVL[:,0],WingShapeAVL[:,1],'r')
    # plt.plot(WingShape[:,0]   ,WingShape[:,1]   ,'b')
    
    # Create AVL File
    if Folder == 0:
        folderName = "Wings\\"+fileName+"\\AVL\\"
    else:
        folderName = Folder+"\\"
    try:
        os.mkdir(folderName)
    except OSError as error:
        print("OVERWRITING EXISTING AVL FILES")
        shutil.rmtree(folderName)
        os.mkdir(folderName)
    
    AVLwingFile = open(folderName+fileName+"Wing.avl", "w")
    
    # Write Header
    AVLwingFile.write(fileName+'\n')
    AVLwingFile.write('#Mach\n')
    AVLwingFile.write('0.0\n')
    AVLwingFile.write('#IYsym   IZsym   Zsym\n')
    AVLwingFile.write(' 0       0       0.0\n')
    AVLwingFile.write('#Sref    Cref    Bref\n')
    np.savetxt(AVLwingFile, [Sref], fmt='%.5f', newline='', header='', footer='    ', comments='', encoding=None)
    np.savetxt(AVLwingFile, [Cref], fmt='%.5f', newline='', header='', footer='   ', comments='', encoding=None)
    np.savetxt(AVLwingFile, [Bref], fmt='%.5f', newline='', header='', footer='\n', comments='', encoding=None)
    AVLwingFile.write('#Xref    Yref    Zref\n')
    AVLwingFile.write(' 0       0       0\n#\n#')
    
    # Write Wing Shape
    NChordwise = nAVL[1]
    Cspace = 1.0
    Nspanwise = nAVL[0]
    Sspace = 1.0
    WingAngle = 0.0
    
    AVLwingFile.write("""
#====================================================================
SURFACE 
Wing 
#Nchordwise  Cspace   Nspanwise   Sspace\n""")
    np.savetxt(AVLwingFile, [NChordwise], fmt='%.0f', newline='', header='', footer='      ', comments='', encoding=None)
    np.savetxt(AVLwingFile, [Cspace], fmt='%.1f', newline='', header='      ', footer='  ', comments='', encoding=None)
    np.savetxt(AVLwingFile, [Nspanwise], fmt='%.0f', newline='', header='    ', footer='    ', comments='', encoding=None)
    np.savetxt(AVLwingFile,    [Sspace], fmt='%.1f', newline='', header='      ', footer='\n', comments='', encoding=None)
    AVLwingFile.write("""#
ANGLE\n""")
    np.savetxt(AVLwingFile,    [WingAngle], fmt='%.1f', newline='', header='', footer='\n', comments='', encoding=None)
    for j in range(nAVL[0]):
        AVLwingFile.write("""#-------------------------------------------------------------
SECTION
#Xle    Yle    Zle     Chord   Ainc  Nspanwise  Sspace\n""")
        np.savetxt(AVLwingFile, [-Bref*xcBV*WingShapeAVL[j][1]], fmt='%.5f', newline='', header='', footer='    ', comments='', encoding=None)
        np.savetxt(AVLwingFile, [Bref*WingShapeAVL[j][0]], fmt='%.5f', newline='', header='', footer='    ', comments='', encoding=None)
        np.savetxt(AVLwingFile, [0], fmt='%.5f', newline='', header='', footer='   ', comments='', encoding=None)
        np.savetxt(AVLwingFile, [Bref*WingShapeAVL[j][1]], fmt='%.5f', newline='', header='', footer='   ', comments='', encoding=None)
        np.savetxt(AVLwingFile, [(180/pi)*WingShapeAVL[j][2]], fmt='%.5f', newline='', header='', footer='   ', comments='', encoding=None)
        np.savetxt(AVLwingFile, [0], fmt='%.5f', newline='', header='', footer='   ', comments='', encoding=None)
        np.savetxt(AVLwingFile, [0], fmt='%.5f', newline='', header='', footer='\n', comments='', encoding=None)
        # add \r\n 
    
    
    AVLwingFile.close()
    
def AVLanalyze(fileName, cL, b, n):
    """
    Inports:
        * filename
            type: string
            description: folder withing "Wings" directory for current wing
        * cL
            type: float
            description: cL at which AVL analyzes the wing
        * b
            type: float
            description: wing span
        * n
            type: int
            description: number for spanwise locations for AVL wing
            
    Returns:
        * AVLdata
            type: n x 13 numpy ndarray
            description:
                rows correspond to spanwise location
                columns correspond to engineering values:
                    col  0: j
                    col  1: Yle
                    col  2: Chord
                    col  3: Area
                    col  4: c cl
                    col  5: ai
                    col  6: cl_norm
                    col  7: cl
                    col  8: cd
                    col  9: cdv
                    col 10: cm_c/4
                    col 11: cm_LE
                    col 12: C.P.x/c
                    
    Suggested feature addition: 
        * ability to read Bref value from input file
    """
    
    import os
    import subprocess as sp
    import numpy as np
    
    AVLexeFolder = "CodeFiles\\Aerodynamics\\Analysis\\"
    AVLwingFolder = "Wings\\"+fileName+"\\AVL\\"
    # AVLanalysisFile = fileName+"AVLanalysisCL"+str(cL)+".txt"
    
    # # Write Commands
    # commandString0 = "LOAD " + AVLwingFolder + fileName + "Wing.avl\n"
    # commandString1 = "oper\n"
    # commandString2 = "A C "+str(cL)+"\n"
    # commandString3 = "X\n"
    # commandString4 = "FS "+AVLanalysisFile+"\n \n"
    # commandString5 = "quit\n"
    AVLstripAnalysisFile = fileName+"AVLstripAnalysis"+str(cL)+".txt"
    AVLtotalAnalysisFile = fileName+"AVLtotalAnalysis"+str(cL)+".txt"
    
    # Write Commands
    commandString0 = "LOAD " + AVLwingFolder + fileName + "Wing.avl\n"
    commandString1 = "oper\n"
    commandString2 = "A C "+str(cL)+"\n"
    commandString3 = "X\n"
    commandString4 = "FS "+AVLstripAnalysisFile+"\n"
    commandString5 = "FT "+AVLtotalAnalysisFile+"\n \n"
    commandString6 = "quit\n"
    
    print("AVL Analysis Started@ cL = "+str(cL))
    
    # Start subprocess to run AVL
    ps = sp.Popen([AVLexeFolder + "avl.exe"],
                  stdin=sp.PIPE,stdout=None,stderr=None)
    # Output commands to AVL
    ps.stdin.write(bytes(commandString0, 'utf-8'))
    ps.stdin.write(bytes(commandString1, 'utf-8'))
    ps.stdin.write(bytes(commandString2, 'utf-8'))
    ps.stdin.write(bytes(commandString3, 'utf-8'))
    ps.stdin.write(bytes(commandString4, 'utf-8'))
    ps.stdin.write(bytes(commandString5, 'utf-8'))
    ps.stdin.write(bytes(commandString6, 'utf-8'))
    # 
    
    # Close subprocess after 60 sec.  This forces output txt file to appear.
    try:
        output, error = ps.communicate(timeout=60)  # Stops trying after 1 min
        print("AVL Analysis FINISHED")
    except sp.TimeoutExpired:
        ps.kill()
        print("XXXXXXXX AVL TIMEOUT ERROR!! XXXXXXXX")
        # output, error = ps.communicate()
        # print(output)
        # print(error)
    
    
    AVLvectors = np.loadtxt(AVLstripAnalysisFile,skiprows = 20,max_rows = n, usecols=(0, 1,2,3,4,5,6,7,8,9,10,11))
    
    # b = max(AVLvectors[:,1]) - min(AVLvectors[:,1])
    
    # fo = open(AVLanalysisFile, "r+")
    # cLfound = 0
    # cDfound = 0
    # Sfound = 0
    # AVLscalars = [0,0,0]
    # for i in range(20):
    #     line = fo.readline()
    #     Sfind = line.find('Surface area')
    #     cLfind = line.find('CLsurf')
    #     cDfind = line.find('CDisurf')
    #     if Sfind != -1:
    #         S    = float(line[Sfind+15:Sfind+26])
    #         cAve = float(line[Sfind+46:Sfind+57])
    #         # b    = S/cAve 
    #         Sfound = 1
    #     if cLfind != -1:
    #         AVLscalars[0] = float(line[cLfind+10:cLfind+19])
    #         cLfound = 1
    #     if cDfind != -1:
    #         AVLscalars[1] = float(line[cDfind+10:cDfind+19])
    #         cDfound = 1
    #     if cLfound == 1 and cDfound == 1 and Sfound == 1:
    #         break
    # fo.close()
    # if Sfound == 0:
    #     print('ERROR: S NOT FOUND !!!')
    # if cLfound == 0:
    #     print('ERROR: cL NOT FOUND !!!')
    # if cDfound == 0:
    #     print('ERROR: cD NOT FOUND !!!')
    # if cLfound and Sfound:
    #     AR = b**2/S
    #     AVLscalars[2] = AVLscalars[0]**2/(np.pi*AVLscalars[1]*AR)
    
    # print("AVL: cL  = ",AVLscalars[0])
    # print("AVL: cDi = ",AVLscalars[1])
    # print("AVL: e   = ",AVLscalars[2])
    # print("AVL: AR  = ",AR)
        
    # # Move analysis text file to current wing's folder
    # os.rename(AVLanalysisFile, AVLwingFolder+AVLanalysisFile)
    
    fo = open(AVLtotalAnalysisFile, "r+")
    cLtotFound = 0
    cLffFound = 0
    cDiindFound = 0
    cDiffFound = 0
    eFound = 0
    for i in range(30):
        line = fo.readline()
        cLtotFind  = line.find('CLtot')
        cLffFind   = line.find('CLff')
        cDiindFind = line.find('CDind')
        cDiffFind  = line.find('CDff')
        eFind      = line.find('  e ')
        
        if cLtotFind != -1:
            CLtot = float(line[cLtotFind+10:cLtotFind+20])
            cLtotFound = 1
        if cLffFind != -1:
            CLff = float(line[cLffFind+10:cLffFind+20])
            cLffFound = 1     
        if cDiindFind != -1:
            cDiind = float(line[cDiindFind+8:cDiindFind+18])
            cDiindFound = 1
        if cDiffFind != -1:
            cDiff = float(line[cDiffFind+8:cDiffFind+18])
            cDiffFound = 1
        if eFind != -1:
            e_span = float(line[eFind+8:eFind+18])
            eFound = 1
        if cLtotFound == 1 and cLffFound == 1 and cDiindFound and cDiffFound and eFound:
            break
    fo.close()
    if cLtotFound == 0:
        print('ERROR: cLtot NOT FOUND !!!')
    if cLffFound == 0:
        print('ERROR: cLff NOT FOUND !!!')
    if cDiindFound == 0:
        print('ERROR: cDind NOT FOUND !!!')
    if cDiffFound == 0:
        print('ERROR: cDiff NOT FOUND !!!')
    if eFound == 0:
        print('ERROR: e NOT FOUND !!!')
    
    AVLscalars = [CLtot,cDiind,e_span]
    # AVLscalars = [CLff,cDiff,e_span]
    print("AVL: cL  = ",AVLscalars[0])
    print("AVL: cDi = ",AVLscalars[1])
    print("AVL: e   = ",AVLscalars[2])
    # print("AVL: AR  = ",AR)
        
    # Move analysis text file to current wing's folder
    os.rename(AVLstripAnalysisFile, AVLwingFolder+AVLstripAnalysisFile)
    os.rename(AVLtotalAnalysisFile, AVLwingFolder+AVLtotalAnalysisFile)
    
    return (AVLvectors,AVLscalars)

def AVLcompare(cL, WingShape,A, fileName, b,nAVL, xcBV, dispPlots=1):
    """
    Inputs:
        * WingShape
            type: n x 3 numpy array
            description: defines wing shape
                rows correspond to spanwise location
                columns correspond to engineering values:
                    col  0: y
                    col  1: c/b
                    col  2: alpha [rad]
        * A
            type:  float or matrix
            description: defines sectional lift curve slope.  If float, then 
                constant value for entire span. If matrix, then spanwise 
                distribution
        * alpha_w
            type: float or array
            description: alpha_w values at which to evaluate wing
        * disp plot
            type: integer (0 or 1)
            description: turns plotting on (1) or off (0)
            
    Returns:
        * vectors
            type:  n x 3 x ncL numpy array
            description: contain the spanwise distributions of 
                (col 0) nonDim Vorticity, (col 1) cl, and (col 2) alpha_i for 
                each cL condition
        * scalars
            type:  3 x ncL numpy array
            description: contain the scalar values of (row 0) cL, (row 1) cDi, 
                and (row 2) span efficiency for each cL condition
    """
    import numpy as np
    n = len(WingShape[:,0])
    AR = n/sum(WingShape[:,1])
    
    AVLwrite(fileName,b,AR,nAVL,xcBV,WingShape)
    ncL = np.size(cL)
    AVLvectors = np.zeros((nAVL[0],12,ncL))
    AVLscalars = np.zeros((3,ncL))
    for i in range(ncL):
        [AVLvectors[:,:,i],AVLscalars[:,i]] = AVLanalyze(fileName, cL[i], b, nAVL[0])
    [LLvectors,LLscalars] = LLanalyzeCL(WingShape,A,cL,0)
    
    if dispPlots:
    
        import matplotlib.pyplot as plt
        lineStyleOrder = ['--',':','-.']
        
        width  = 6.5
        height = 6.5
        plt.figure(figsize=(width, height))
        plt.subplot(311)
        plt.grid()
        plt.ylabel(r'$c c_\ell/c_\mathrm{ref}$', fontsize=20)
        # plt.axis('equal')
        maxcCL = 0
        for i in range(ncL):
            plt.plot(AVLvectors[:,1,ncL-1-i]/b, AVLvectors[:,4,ncL-1-i]*AR/b, 'r'+lineStyleOrder[i%3],
                      # label = "AVL:$c_L$=%3.2f; $c_{Di}$=%5.4f"%(AVLscalars[0,i],AVLscalars[1,i]))
                       label = "AVL:cL=%3.2f; cDi=%5.4f; e=%5.4f"%(AVLscalars[0,ncL-1-i],AVLscalars[1,ncL-1-i],AVLscalars[2,ncL-1-i]))
            cCLcref = LLvectors[:,1,ncL-1-i]*WingShape[:,1]*AR
            plt.plot(WingShape[:,0], cCLcref, 'k'+lineStyleOrder[i%3],
                       label = " LL:cL=%3.2f; cDi=%5.4f; e=%5.4f"%(LLscalars[0,ncL-1-i],LLscalars[1,ncL-1-i],LLscalars[2,ncL-1-i]))
            maxcCL = max(maxcCL,max(cCLcref))
        plt.xticks(np.array((-0.5,-0.25,0,0.25,0.5)), (), fontsize=18)
        plt.yticks(fontsize=18)
        plt.xlim([-0.5,0.5])
        plt.ylim(top=1.1*maxcCL)
        plt.legend(fontsize=18, loc='lower center',bbox_to_anchor=(0.5, 1.02),ncol=1);
        
        plt.subplot(312)
        plt.grid()
        plt.ylabel(r'$c_\ell$', fontsize=20)
        # plt.axis('equal')
        for i in range(ncL):
            plt.plot(AVLvectors[:,1,ncL-1-i]/b, AVLvectors[:,7,ncL-1-i], 'r'+lineStyleOrder[i%3])
            plt.plot(WingShape[:,0], LLvectors[:,1,ncL-1-i],'k'+lineStyleOrder[i%3])
        plt.xticks(np.array((-0.5,-0.25,0,0.25,0.5)), (), fontsize=18)
        plt.yticks(fontsize=18)
        plt.xlim([-0.5,0.5])
        # plt.xlabel(r'$y/b$', fontsize=18)
        
        plt.subplot(313)
        plt.grid()
        plt.ylabel(r'$[\alpha_i]_\mathrm{ff}, deg$', fontsize=18)
        # plt.axis('equal')
        minAlphai = 0
        for i in range(ncL):
            # plt.plot(AVLvectors[:,1,i]/b, AVLvectors[:,5,i]*180/np.pi, 'r'+lineStyleOrder[i%3])
            plt.plot(AVLvectors[:,1,ncL-1-i]/b, AVLvectors[:,5,ncL-1-i]*180/np.pi, 'r'+lineStyleOrder[i%3])
            plt.plot(     WingShape[:,0],  2*LLvectors[:,2,ncL-1-i]*180/np.pi, 'k'+lineStyleOrder[i%3])
            minAlphai = min(min(AVLvectors[:,5,ncL-1-i]),min(LLvectors[:,2,ncL-1-i]))
        plt.xticks(np.array((-0.5,-0.25,0,0.25,0.5)), fontsize=18)
        # plt.yticks(np.array((-5,-4,-3,-2,-1,0,1,2,3,4,5)))
        plt.yticks(fontsize=18)
        plt.xlim([-0.5,0.5])
        # if minAlphai >= 0:
        #     plt.ylim(bottom=0)
        plt.xlabel(r'$y/b$', fontsize=20)
    
    return(LLvectors,LLscalars)

def AVLanalyzeAlpha(fileName, alpha, b, n):
    """
    Inports:
        * filename
            type: string
            description: folder withing "Wings" directory for current wing
        * alpha
            type: float
            description: angle of attack in degrees at which AVL analyzes 
                         the wing
        * b
            type: float
            description: wing span
        * n
            type: int
            description: number for spanwise locations for AVL wing
            
    Returns:
        * AVLvectors
            type: n x 13 numpy ndarray
            description:
                rows correspond to spanwise location
                columns correspond to engineering values:
                    col  0: j
                    col  1: Yle
                    col  2: Chord
                    col  3: Area
                    col  4: c cl
                    col  5: ai
                    col  6: cl_norm
                    col  7: cl
                    col  8: cd
                    col  9: cdv
                    col 10: cm_c/4
                    col 11: cm_LE
                    col 12: C.P.x/c
        * AVLscalars
            type: 3 element list
            description:  Scalars of wing [cL, cDi, e]

                    
    Suggested feature addition: 
        * ability to read Bref value from input file
    """
    
    import os
    import subprocess as sp
    import numpy as np
    
    AVLexeFolder = "CodeFiles\\Aerodynamics\\Analysis\\"
    AVLwingFolder = "TipStudy\\"
    AVLstripAnalysisFile = fileName+"AVLstripAnalysisAlpha"+str(alpha)+".txt"
    AVLtotalAnalysisFile = fileName+"AVLtotalAnalysisAlpha"+str(alpha)+".txt"
    
    # Write Commands
    commandString0 = "LOAD " + AVLwingFolder + fileName + "Wing.avl\n"
    commandString1 = "oper\n"
    commandString2 = "A A "+str(alpha)+"\n"
    commandString3 = "X\n"
    commandString4 = "FS "+AVLstripAnalysisFile+"\n"
    commandString5 = "FT "+AVLtotalAnalysisFile+"\n \n"
    commandString6 = "quit\n"
    # commandString6 = "exit\n"
    
    print("AVL Analysis Started@ alpha = "+str(alpha))
    
    # Start subprocess to run AVL
    ps = sp.Popen([AVLexeFolder + "avl.exe"],
                  stdin=sp.PIPE,stdout=None,stderr=None)
    # Output commands to AVL
    ps.stdin.write(bytes(commandString0, 'utf-8'))
    ps.stdin.write(bytes(commandString1, 'utf-8'))
    ps.stdin.write(bytes(commandString2, 'utf-8'))
    ps.stdin.write(bytes(commandString3, 'utf-8'))
    ps.stdin.write(bytes(commandString4, 'utf-8'))
    ps.stdin.write(bytes(commandString5, 'utf-8'))
    ps.stdin.write(bytes(commandString6, 'utf-8'))
    # ps.stdin.write(bytes(commandString6, 'utf-8'))
    # 
    
    # Close subprocess after 60 sec.  This forces output txt file to appear.
    try:
        output, error = ps.communicate(timeout=60)  # Stops trying after 1 min
        print("AVL Analysis FINISHED")
    except sp.TimeoutExpired:
        ps.kill()
        print("XXXXXXXX AVL TIMEOUT ERROR!! XXXXXXXX")
        # output, error = ps.communicate()
        # print(output)
        # print(error)
    
    
    AVLvectors = np.loadtxt(AVLstripAnalysisFile,skiprows = 20,max_rows = n)
    
    # b = max(AVLvectors[:,1]) - min(AVLvectors[:,1])
    
    # fo = open(AVLstripAnalysisFile, "r+")
    # cLfound = 0
    # cDfound = 0
    # Sfound = 0
    # AVLscalars = [0,0,0]
    # for i in range(20):
    #     line = fo.readline()
    #     Sfind = line.find('Surface area')
    #     cLfind = line.find('CLsurf')
    #     cDfind = line.find('CDisurf')
    #     if Sfind != -1:
    #         S    = float(line[Sfind+15:Sfind+26])
    #         cAve = float(line[Sfind+46:Sfind+57])
    #         # b    = S/cAve 
    #         Sfound = 1
    #     if cLfind != -1:
    #         AVLscalars[0] = float(line[cLfind+10:cLfind+19])
    #         cLfound = 1
    #     if cDfind != -1:
    #         AVLscalars[1] = float(line[cDfind+10:cDfind+19])
    #         cDfound = 1
    #     if cLfound == 1 and cDfound == 1 and Sfound == 1:
    #         break
    # fo.close()
    # if Sfound == 0:
    #     print('ERROR: S NOT FOUND !!!')
    # if cLfound == 0:
    #     print('ERROR: cL NOT FOUND !!!')
    # if cDfound == 0:
    #     print('ERROR: cD NOT FOUND !!!')
    # if cLfound and Sfound:
    #     AR = b**2/S
    #     AVLscalars[2] = AVLscalars[0]**2/(np.pi*AVLscalars[1]*AR)
    #     print("calculated e!!!!!") 
    
    
    fo = open(AVLtotalAnalysisFile, "r+")
    cLtotFound = 0
    cLffFound = 0
    cDiindFound = 0
    cDiffFound = 0
    eFound = 0
    for i in range(30):
        line = fo.readline()
        cLtotFind  = line.find('CLtot')
        cLffFind   = line.find('CLff')
        cDiindFind = line.find('CDind')
        cDiffFind  = line.find('CDff')
        eFind      = line.find('  e ')
        
        if cLtotFind != -1:
            CLtot = float(line[cLtotFind+10:cLtotFind+20])
            cLtotFound = 1
        if cLffFind != -1:
            CLff = float(line[cLffFind+10:cLffFind+20])
            cLffFound = 1     
        if cDiindFind != -1:
            cDiind = float(line[cDiindFind+8:cDiindFind+18])
            cDiindFound = 1
        if cDiffFind != -1:
            cDiff = float(line[cDiffFind+8:cDiffFind+18])
            cDiffFound = 1
        if eFind != -1:
            e_span = float(line[eFind+8:eFind+18])
            eFound = 1
        if cLtotFound == 1 and cLffFound == 1 and cDiindFound and cDiffFound and eFound:
            break
    fo.close()
    if cLtotFound == 0:
        print('ERROR: cLtot NOT FOUND !!!')
    if cLffFound == 0:
        print('ERROR: cLff NOT FOUND !!!')
    if cDiindFound == 0:
        print('ERROR: cDind NOT FOUND !!!')
    if cDiffFound == 0:
        print('ERROR: cDiff NOT FOUND !!!')
    if eFound == 0:
        print('ERROR: e NOT FOUND !!!')
    
    # AVLscalars = [CLtot,cDiind,e_span]
    AVLscalars = [CLff,cDiff,e_span]
    print("AVL: cL  = ",AVLscalars[0])
    print("AVL: cDi = ",AVLscalars[1])
    print("AVL: e   = ",AVLscalars[2])
    # print("AVL: AR  = ",AR)
        
    # Move analysis text file to current wing's folder
    os.rename(AVLstripAnalysisFile, AVLwingFolder+AVLstripAnalysisFile)
    os.rename(AVLtotalAnalysisFile, AVLwingFolder+AVLtotalAnalysisFile)
    
    return (AVLvectors,AVLscalars)