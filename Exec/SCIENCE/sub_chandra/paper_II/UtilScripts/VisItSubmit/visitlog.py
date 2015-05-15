# Visit 2.8.0 log file
ScriptVersion = "2.8.0"
if ScriptVersion != Version():
    print "This script is for VisIt %s. It may not work with version %s" % (ScriptVersion, Version())
ShowAllWindows()
OpenDatabase("/u/sciteam/ajacobs/scratch/sub_chandra/11030-108-185-4lev/plotfiles/11030-108-185-4lev_plt08568/Header", 0)
# The UpdateDBPluginInfo RPC is not supported in the VisIt module so it will not be logged.
AddPlot("Contour", "X_lp_C12_rp_", 1, 1)
DrawPlots()
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 1
SaveWindowAtts.outputDirectory = "."
SaveWindowAtts.fileName = "visit"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.JPEG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
SaveWindowAtts.width = 1024
SaveWindowAtts.height = 1024
SaveWindowAtts.screenCapture = 0
SaveWindowAtts.saveTiled = 0
SaveWindowAtts.quality = 80
SaveWindowAtts.progressive = 0
SaveWindowAtts.binary = 0
SaveWindowAtts.stereo = 0
SaveWindowAtts.compression = SaveWindowAtts.PackBits  # None, PackBits, Jpeg, Deflate
SaveWindowAtts.forceMerge = 0
SaveWindowAtts.resConstraint = SaveWindowAtts.ScreenProportions  # NoConstraint, EqualWidthHeight, ScreenProportions
SaveWindowAtts.advancedMultiWindowSave = 0
SetSaveWindowAttributes(SaveWindowAtts)
SaveWindow()
ClearCacheForAllEngines()
SetActivePlots(0)
DeleteActivePlots()
OpenDatabase("/u/sciteam/ajacobs/scratch/sub_chandra/11030-108-185-4lev/plotfiles/11030-108-185-4lev_plt08568/Header")
ClearCacheForAllEngines()
CloseComputeEngine("nid27638", "")
