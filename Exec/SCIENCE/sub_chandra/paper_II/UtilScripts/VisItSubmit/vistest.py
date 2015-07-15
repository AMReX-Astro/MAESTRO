#This is intended to be sourced by a batch script submitted via qssub.  
#The VisIt python library will be imported.  This is where any functions
#that are called yet aren't imported come from.

#Usage: visit -cli -nowin -s vistest.py

### Imports
from sys import exit

### Execution
#Configure window properties
s = SaveWindowAttributes()
s.format = s.JPEG
s.width, s.height = 1024,1024
SetSaveWindowAttributes(s)

#Open the file/database
db = "/u/sciteam/ajacobs/scratch/sub_chandra/11030-108-185-4lev/plotfiles/11030-108-185-4lev_plt08568/Header"
OpenDatabase(db)

#Configure annotation

#Add plots
AddPlot('Contour', 'X_lp_C12_rp_')

#Configure plots

#Draw and save
DrawPlots()
SaveWindow()

#Clean up, finalize
ClearCacheForAllEngines()
DeleteAllPlots()
CloseDatabase(db)
ClearCacheForAllEngines()
CloseComputeEngine()
sys.exit() #This is required when sourcing a script into VisIt's CLI
