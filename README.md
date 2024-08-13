# Read Me: G-Wing

Software Name:  G-Wing
Author:         Justin Valenti

## Summary:
G-Wing is a tool developed to the streamline the design process of wings built with Fused Filament Fabrication.  G-Wing was developed as part of Justin's doctoral work.

## Running G-Wing:
Note: G-Wing was developed in the Python 3 language, using the Spyder IDE.

1. Open "Gwing/GWingMAIN.py"
2. In the text of GWingMAIN, uncomment the desired wing 
    * Preloaded wings are commented out in a listed immediately following header comments
    * One wing will be uncommented, so G-Wing will execute without error with no editing.
3. Execute GWingMAIN.py
4. View plot outputs to view the wing being designed.  See "Edditing Existing Wings" to learn how to output G-Code

## Editing Existing Wings (And outputing G-Code)
1. Open "GWing/Wings/<InsertWingName>/Settings.pi"
2. Edited desired parameter
    * The attempt was made to thoroughly comment the settings file
    * A typical parameter to edit is "GCodeWrite", this comes default as 0.  Change this parameter to 1 for G-Wing to generate G-Code.  Once executed, the resulting G-Code will be placed in the "GWing/Wings/<InsertWingName>/GCode" folder.

## Creating New Wings
1. In the "Wings" folder, copy any existing wing folder and paste the copy as the new name, for now lets call this new wing "NewWing".
2. Optional: Delete the "AVL" and "GCode" folders---new folders will be generated when the wing is generated the first time.
3. In "GWing/Wings/NewWing/Settings.py", update the fileName string to the new name
4. In GWingMAIN, add "from Wings.NewWing.Settings import *" to the bottom of the list.  Leave the line uncommented to generate the new wing, comment the line out during execution to not generate this wing.

## Included OpenSource Software

G-Wing includes two opensource software:
1. AVL (Athena Vortex Lattice) is an open-source aerodynamic analysis code available here: 
http://web.mit.edu/drela/Public/web/avl/
2. OpenSCAD is an open source CAD software available here:
https://openscad.org/downloads.html