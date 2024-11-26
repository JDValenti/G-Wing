# Collection of Python functions to directly manipulate GCode
# Special thanks to Natalie Nash for help in developing Mirror


def Mirror(FileName):
    import re
    
    # Read in GCode, measure length
    oldFile = open(FileName+".gcode","r")
    ExistingGC = oldFile.readlines(  )
    oldFile.close()
    LenGC = len(ExistingGC)
    
    # Identify Body of GCode, Measure Max & Min X Coords
    inBody = 0
    xMin = None
    xMax = None
    for i in range(LenGC):
        if inBody:
            if ExistingGC[i] == ";Begin Custom End GCode\n":
                inBody = 0
                break
            xMatch = re.search(r' X(.*?) ', ExistingGC[i])
            if xMatch:
                x = float(xMatch.group(1))
                if xMin == None:
                    xMin = x
                    xMax = x
                else:
                    if x < xMin:
                        xMin = x
                    elif x > xMax:
                        xMax = x

        else:
            if ExistingGC[i] == ";Finished Custom Start GCode\n":
                inBody = 1
    
    # Begin second sweep to modify X coordinates
    inBody = 0
    newFile = open(FileName + "_M.gcode", "w") # Open new file
    for i in range(LenGC): # Iterate through each line of GCode
        if inBody:
            if ExistingGC[i] == ";Begin Custom End GCode\n": # Check if at footer
                inBody = 0
                continue
            newString = ""
            xMatch = re.search(r' X(.*?) ', ExistingGC[i]) # Check if 
            if xMatch:
                xStr = xMatch.group(0)[2:-1]
                xPrcsn = len(xStr.split('.')[1])
                x = float(xMatch.group(1))
                xNew = round(xMin + xMax - x, xPrcsn)
                xNewStr = str(round(xNew, xPrcsn))
                xNewPrcsn = len(xNewStr.split('.')[1])
                dPrcsn = xPrcsn - xNewPrcsn
                for j in range(dPrcsn):
                    xNewStr = xNewStr + "0"
                oldString = ExistingGC[i]
                newString = oldString.replace("X"+xStr, "X"+xNewStr)
            else:
                newString = ExistingGC[i]
            newFile.write(newString)

        else:
            newFile.write(ExistingGC[i])
            # Check if finished with header
            if ExistingGC[i] == ";Finished Custom Start GCode\n": 
                inBody = 1
                

    newFile.close()
    print("Mirrored!")
    return "done"

def Center(FileName,xCenter,yCenter):
    import re
    
    oldFile = open(FileName+".gcode","r")
    ExistingGC = oldFile.readlines(  )
    oldFile.close()
    LenGC = len(ExistingGC)
    
    # Identify Body of GCode, Measure Max & Min X Coords
    inBody = 0
    xMin = None
    xMax = None
    yMin = None
    yMax = None
    for i in range(LenGC):
        if inBody:
            if ExistingGC[i] == ";Begin Custom End GCode\n":
                inBody = 0
                break
            
            # Find min/max x coordinate
            xMatch = re.search(r' X(.*?) ', ExistingGC[i])
            if xMatch:
                x = float(xMatch.group(1))

                if xMin == None:
                    xMin = x
                    xMax = x
                else:
                    if x < xMin:
                        xMin = x
                    elif x > xMax:
                        xMax = x
            # Find min/max y coordinate
            yMatch = re.search(r' Y(.*?) ', ExistingGC[i])
            if yMatch:
                y = float(yMatch.group(1))

                if yMin == None:
                    yMin = y
                    yMax = y
                else:
                    if y < yMin:
                        yMin = y
                    elif y > yMax:
                        yMax = y
        else:
            if ExistingGC[i] == ";Finished Custom Start GCode\n":
                inBody = 1
    xBarOld = 0.5*(xMin + xMax)
    yBarOld = 0.5*(yMin + yMax)
    
    xTranslate = xCenter - xBarOld
    yTranslate = yCenter - yBarOld

    inBody = 0
    newFile = open(FileName+"_C.gcode", "w")
    for i in range(LenGC):
        if inBody:
            
            if ExistingGC[i] == ";Begin Custom End GCode\n":
                inBody = 0
                continue
            
            newString = ""
            xMatch = re.search(r' X(.*?) ', ExistingGC[i])
            if xMatch:
                xStr = xMatch.group(0)[2:-1]
                xPrcsn = len(xStr.split('.')[1])
                x = float(xMatch.group(1))
                xNew = round(x + xTranslate, xPrcsn)
                xNewStr = str(round(xNew, xPrcsn))
                xNewPrcsn = len(xNewStr.split('.')[1])
                dPrcsn = xPrcsn - xNewPrcsn
                for j in range(dPrcsn):
                    xNewStr = xNewStr + "0"
                oldString = ExistingGC[i]
                newString = oldString.replace("X"+xStr, "X"+xNewStr)
            else:
                newString = ExistingGC[i]
                
            yMatch = re.search(r' Y(.*?) ', newString)
            if yMatch:
                yStr = yMatch.group(0)[2:-1]
                yPrcsn = len(yStr.split('.')[1])
                y = float(yMatch.group(1))
                yNew = round(y + yTranslate, yPrcsn)
                yNewStr = str(round(yNew, yPrcsn))
                yNewPrcsn = len(yNewStr.split('.')[1])
                dPrcsn = yPrcsn - yNewPrcsn
                for j in range(dPrcsn):
                    yNewStr = yNewStr + "0"
                oldString = newString
                newString = oldString.replace("Y"+xStr, "Y"+xNewStr)
            newFile.write(newString)

        else:
            newFile.write(ExistingGC[i])
            if ExistingGC[i] == ";Finished Custom Start GCode\n":
                inBody = 1
                

    newFile.close()
    print("Centered!")
    return "done"

def Translate(FileName,xTranslate,yTranslate):
    import re
    
    oldFile = open(FileName+".gcode","r")
    ExistingGC = oldFile.readlines(  )
    oldFile.close()
    LenGC = len(ExistingGC)
    
    inBody = 0
    newFile = open(FileName+"_T.gcode", "w")
    for i in range(LenGC):
        if inBody:
            
            if ExistingGC[i] == ";Begin Custom End GCode\n":
                inBody = 0
                continue
            
            newString = ""
            xMatch = re.search(r' X(.*?) ', ExistingGC[i])
            if xMatch:
                xStr = xMatch.group(0)[2:-1]
                xPrcsn = len(xStr.split('.')[1])
                x = float(xMatch.group(1))
                xNew = round(x + xTranslate, xPrcsn)
                xNewStr = str(round(xNew, xPrcsn))
                xNewPrcsn = len(xNewStr.split('.')[1])
                dPrcsn = xPrcsn - xNewPrcsn
                for j in range(dPrcsn):
                    xNewStr = xNewStr + "0"
                oldString = ExistingGC[i]
                newString = oldString.replace("X"+xStr, "X"+xNewStr)
            else:
                newString = ExistingGC[i]
                
            yMatch = re.search(r' Y(.*?) ', newString)
            if yMatch:
                y = float(yMatch.group(1))
                digitMatch = re.search(r'\.(.*)', newString)
                digitCount = 0
                if digitMatch:
                    digitCount = len(digitMatch.group(1))
                newY = y + yTranslate
                oldString = newString
                newString = oldString.replace(str(y), str(round(newY, digitCount)))
            newFile.write(newString)

        else:
            newFile.write(ExistingGC[i])
            if ExistingGC[i] == ";Finished Custom Start GCode\n":
                inBody = 1
                

    newFile.close()
    return "done"

def Rotate(FileName,theta):
    import numpy as np
    import re
    
    oldFile = open(FileName+".gcode","r")
    ExistingGC = oldFile.readlines(  )
    oldFile.close()
    LenGC = len(ExistingGC)
    
    # Identify Body of GCode, Measure Max & Min X Coords
    inBody = 0
    xMin = None
    xMax = None
    yMin = None
    yMax = None
    for i in range(LenGC):
        if inBody:
            if ExistingGC[i] == ";Begin Custom End GCode\n":
                inBody = 0
                break
            
            # Find min/max x coordinate
            xMatch = re.search(r' X(.*?) ', ExistingGC[i])
            if xMatch:
                x = float(xMatch.group(1))

                if xMin == None:
                    xMin = x
                    xMax = x
                else:
                    if x < xMin:
                        xMin = x
                    elif x > xMax:
                        xMax = x
            # Find min/max y coordinate
            yMatch = re.search(r' Y(.*?) ', ExistingGC[i])
            if yMatch:
                y = float(yMatch.group(1))

                if yMin == None:
                    yMin = y
                    yMax = y
                else:
                    if y < yMin:
                        yMin = y
                    elif y > yMax:
                        yMax = y
        else:
            if ExistingGC[i] == ";Finished Custom Start GCode\n":
                inBody = 1
    xBar = 0.5*(xMin + xMax)
    yBar = 0.5*(yMin + yMax)
    
    inBody = 0
    newFile = open(FileName + "_R.gcode", "w")
    for i in range(LenGC):
        if inBody:
            
            if ExistingGC[i] == ";Begin Custom End GCode\n":
                inBody = 0
                continue
            
            newString = ""
            xMatch = re.search(r' X(.*?) ', ExistingGC[i])
            if xMatch:
                xStr = xMatch.group(0)[2:-1]
                xPrcsn = len(xStr.split('.')[1])
                x = float(xMatch.group(1))
                xC = x - xBar
                
            yMatch = re.search(r' Y(.*?) ', ExistingGC[i])
            if yMatch:
                yStr    = yMatch.group(0)[2:-1]
                yPrcsn  = len(yStr.split('.')[1])
                y       = float(yMatch.group(1))
                yC      = y - yBar
                
            if xMatch:
                xR = xC*np.cos((np.pi/180)*theta) - yC*np.sin((np.pi/180)*theta)
                xNew = round(xR + xBar, xPrcsn)
                xNewStr = str(round(xNew, xPrcsn))
                xNewPrcsn = len(xNewStr.split('.')[1])
                dPrcsn = xPrcsn - xNewPrcsn
                for j in range(dPrcsn):
                    xNewStr = xNewStr + "0"
                oldString = ExistingGC[i]
                newString = oldString.replace("X"+xStr, "X"+xNewStr)
            else:
                newString = ExistingGC[i]
                
            if yMatch:
                yR = xC*np.sin((np.pi/180)*theta) + yC*np.cos((np.pi/180)*theta)
                yNew = yR + yBar
                yNewStr = str(round(yNew, yPrcsn))
                yNewPrcsn = len(yNewStr.split('.')[1])
                dPrcsn = yPrcsn - yNewPrcsn
                for j in range(dPrcsn):
                    yNewStr = yNewStr + "0"
                oldString = newString
                newString = oldString.replace("Y"+yStr, "Y"+yNewStr)
            
            newFile.write(newString)

        else:
            newFile.write(ExistingGC[i])
            if ExistingGC[i] == ";Finished Custom Start GCode\n":
                inBody = 1
    
    print("Rotated!")
                

    newFile.close()
    return "done"

def Platter(FileNames,bedDims):
    import re
    numFiles = len(FileNames)
    for j in range(numFiles):
        Center(FileName,xCenter,yCenter)
        Translate(FileName,xTranslate,yTranslate)
    
                

    newFile.close()
    print("Centered!")
    return "done"