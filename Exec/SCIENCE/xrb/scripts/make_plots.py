import ConfigParser
import sys
import getopt
import string

globalParams = {}
globalRuns = []
globalNumPlots = 0

def isStrAnInt(string):
    try: int(string)
    except ValueError: return 0
    else: return 1


def isStrAFloat(string):
    try: float(string)
    except ValueError: return 0
    else: return 1


def isBlank(key):
    """
    returns true if dict does not conatain key OR dict does have key but its
    value is an empty string"""
    global globalParams

    try: globalParams[key]
    except KeyError: return 1
    else:
        if globalParams[key] == "": return 1
        return 0


def getParam(key):
    """
    a little more clean than cumbersome dictionary-key based retrieval
    """
    global globalParams

    if globalParams.has_key(key): return globalParams[key]
    else:
        print "ERROR in getParam: %s is not a key in globalParams" % key
        sys.exit()


def setParam(key, value):
    """
    a little more clean than cumbersome dictionary-key based assignment
    """
    global globalParams
    
    globalParams[key] = value


def getRunInfo(file):
    """
    returns a dictionary with keys of the form <run>.<option> and the 
    corresponding values from the file.  
    if some options are not set, we set some defaults 
    """
    global globalParams
    global globalRuns
    global globalNumPlots

    try: f = open(file,'r')
    except IOError:
        print "ERROR: parameter file %s does not exist" % file
        sys.exit()
    else:
        f.close()

    cp = ConfigParser.ConfigParser()
    cp.optionxform = str
    cp.read(file)

    globalRuns = cp.sections()

    for run in cp.sections():

        for opt in cp.options(run):
        
            value = cp.get(run, opt)

            # save the proper datatype
            if (isStrAnInt(value)):
                setParam(run + '.' + opt, int(value))
            elif (isStrAFloat(value)):
                setParam(run + '.' + opt, float(value))
            else:
                setParam(run + '.' + opt, value)

        # make sure we have names for each run
        if (isBlank(run + '.baseName')):
            setParam(run + '.baseName', run)

        startNum = getParam(run + '.startNum')

        # if endNum or stride is not given, we only plot startNum
        if ( isBlank(run + '.endNum') or isBlank(run + '.stride') ):
            setParam(run + '.endNum', startNum)
            setParam(run + '.stride', 1)

        endNum   = getParam(run + '.endNum'  )
        stride   = getParam(run + '.stride'  )
            
        # calculate the total number of plots; this counts the end file as well
        for i in range(startNum, endNum + stride, stride):
            globalNumPlots += 1


def make_plots(argv):

    global globalRuns
    global globalParams

    usage = """
    ./make_plots.py [-p|--plot_these_files pltfile_list, 
                     -v|--vars var_list]
          files_to_plot.ini

    arguments:

      files_to_plot.ini
            This is the input file that defines which pltfiles are to
            be plotted if the -p|--plot_these_files option is NOT set.
            Default behaviour is to read this .ini file.  It has the format

               [xrb_2d_5cm]
               dataDir  = < relative path (from working directory) to the 
                            directory containing the plotfiles >
               startNum = < # of pltfile to start plotting >
               endNum   = < # of pltfile to end plotting >
               stride   = < how often between startNum and endNum to generate
                            a plot >
               vars     = < which vars that are to be plotted >
               baseName = < name to label the plot images >

            Here, [xrb_2d_5cm] is the name to describe the particular run.  
            There can be many such runs to be plotted, each with their own 
            unique name, following the format of [xrb_2d_5cm] above.

            For a particular run, (e.g. [xrb_2d_5cm]),
               
               dataDir is the relative path to the directory that contains
               the various pltfiles; e.g. for the pltfile 
               ./092908/2d_full_domain_10cm/plt00000/, dataDir would be 
               ./092908/2d_full_domain_10cm

               startNum, endNum, and stride give the range of which pltfiles
               are to be plotted along with how many plots are to be generated.
               For example, suppose dataDir contains the following pltfiles:
                   plt00000/
                   plt10000/
                   plt20000/
                   plt30000/
                   plt40000/
               To generate plots for all of pltfiles one would specify 
                   startNum = 0
                   endNum   = 40000
                   stride   = 10000
               To plot every other file one would specify
                   startNum = 0
                   endNum   = 40000
                   stride   = 20000
               To plot a single file, one needs only specify startNum, i.e.
               endNum defaults to startNum.

               vars is a tuple containing the list of vars to be plotted.
               Default is vars = (\"X_He4_\",) which uses the syntax from VisIt
               for the helium mass fraction.

               baseName is the name used to label the plot images; e.g. 
               plotting the density in pltfile plt03200 with 
               baseName = \"2d_small\" would result in the image file
               \"2d_small_density_03200.png\"  Default is to set baseName 
               equal to the run name; e.g. for the [xrb_2d_5cm] run, we would
               have \"xrb_2d_5cm_density_03200.png\"

    
    options:

      -p|--plot_these_files pltfile_list (currently disabled)
        Instead of using a files_to_plot.ini file, generate plots from the
        pltfiles in pltfile_list.  Images saved in this fashion will have
        a baseName equal to the pltfileName in pltfile_list.  The default
        variable plotted will be Helium mass fraction, unless the -v|--vars
        option is set; see below.

      -v|--vars var_list (currently disabled)
        If this option is set and the -p option is set, then plots of each 
        var in var_list will be generated for each pltfile in pltfile_list.
        Images generated in this fashion will have a baseName which contains
        the pltfileName and the varName.
        If this option is set without the -p option, then var_list will 
        override the vars variable for each run listed in files_to_plot.ini.

    """


    if len(argv) == 0:
        print usage
        sys.exit(2)

    try: 
        opts, next = getopt.getopt(argv, "p:v:", 
                                   ["plot_these_files=",
                                    "vars="])

    except getopt.GetoptError:
        print "Invalid calling sequence"
        print usage
        sys.exit(2)

    plot_filelist = 0
#    plot_vars = "X_He4_"
    plot_vars = "Machnumber"

    for o, a in opts:
        
        if o == "-p" or o == "--plot_these_files":
            print "this feature currently disabled; please use a .ini file"
            print "wanting to plot these files:", a
            sys.exit(2)

        if o == "-v" or o == "--vars":
            print "this feature currently disabled; please use a .ini file"
            print "wanting to plot these vars:", a
            sys.exit(2)

    try:
        runFile = next[0]

    except IndexError:
        print "ERROR: a run file was not specified"
        print usage
        sys.exit(2)

    # -------------------------------------------
    # get the run information and setup defaults
    # -------------------------------------------
    getRunInfo(runFile)

    # ----------------------------------------------------------------------
    # create a list of databases (i.e. Header files) to plot along with the
    # names of the corresponding image filename
    # ----------------------------------------------------------------------
    db_list = []
    image_filenames = []

    for run in globalRuns:

        dataDir = getParam(run + '.dataDir' )
        start   = getParam(run + '.startNum')
        end     = getParam(run + '.endNum'  )
        stride  = getParam(run + '.stride'  )

        # if endNum - startNum is not an integer multiple of stride, then we
        # can't use range(start, end+stride, stride) to capture the end 
        # plotfile so we will add this in separateley
        for fileNum in range(start, end, stride):
            
            db_list.append(dataDir + 'plt' + str(fileNum).rjust(5,'0') \
                               + '/Header')
            image_filenames.append(run + '_' + str(fileNum).rjust(5,'0'))

        # add last plotfile
        db_list.append(dataDir + 'plt' + str(end).rjust(5,'0') + '/Header')
        image_filenames.append(run + '_' + str(end).rjust(5,'0'))
            

    # ------------------------------------
    # set up the plotter and its settings
    # ------------------------------------
    import visit

    # launch visit with no visualization windows shown
    visit.AddArgument("-nowin")
    visit.Launch()

    # load the first file as a template for the others
    visit.OpenDatabase(db_list[0])
    visit.AddPlot("Pseudocolor",plot_vars)

    # get objects to various plot attributes
    p = visit.PseudocolorAttributes()
    a = visit.AnnotationAttributes()
    s = visit.SaveWindowAttributes()

    # change the colortable, range, and scaling
    p.colorTableName="caleblack"
    p.minFlag, p.min = 1, 1.0e-10
    p.maxFlag, p.max = 1, 1.0e0
    p.scaling=p.Log
        
    # add some unit labels to the axes
    a.axes2D.xAxis.title.userUnits, a.axes2D.xAxis.title.units = 1, "cm"
    a.axes2D.yAxis.title.userUnits, a.axes2D.yAxis.title.units = 1, "cm"
    # turn off some unwanted database info
    a.userInfoFlag, a.databaseInfoFlag = 0,0

    # set save format to be png and turn off VisIt's annoying habit of 
    # appending numbers to the end of filenames
    s.format = s.PNG
    s.quality = 100
    s.family = 0
    s.fileName = image_filenames[0]
    s.width, s.height = 600, 600

    # apply the appearance changes
    visit.SetPlotOptions(p)
    visit.SetAnnotationAttributes(a)
    visit.SetSaveWindowAttributes(s)

    # render the plot
    visit.DrawPlots()
        
    # save the image
    name = visit.SaveWindow()

    # open the remaining files 
    for plot in range(1,globalNumPlots):
        # close the old database and add a new one with the same plot
        # attributes as above
        visit.ReplaceDatabase(db_list[plot])

        # change the save name
        s.fileName = image_filenames[plot]
        visit.SetSaveWindowAttributes(s)

        # render the new plot
        visit.DrawPlots()
        
        # save the image
        name = visit.SaveWindow()


if __name__ == "__main__":
    # drop the program name from argv
    make_plots(sys.argv[1:])


