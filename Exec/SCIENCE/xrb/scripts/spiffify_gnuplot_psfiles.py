#! /usr/bin/python

import sys
import re
import getopt


usage = """
 spiffify_gnuplot_psfiles.py [-012345678 color -t thickness -s] <file1 ...>

   Description:
     This routine applies changes to PostScript files created by gnuplot with
     the "set terminal postscript color enhanced" option.  

     -012345678 r:g:b          Specify the color for the particular linetype.
                               Note that linetypes are numbered starting from 
                               zero.  Currently, the color should be specified
                               in an r:g:b format where r,g,b range from 0 to 
                               255 and correspond to the red, green, and blue
                               primary colors.  If not specified, the defaults
                               from gnuplot will be used.

     -t thickness              Specify the thickness of the plotted lines.
                               gnuplot scales the unit for this to be 
                               1/720th of an inch, or 1/10th of a pt.  Default
                               value is 15.

     -s                        Specify that all lines should be of solid 
                               linetype.
                               
"""

userLineWidthRE = re.compile('/userlinewidth')
PLRE = re.compile('/PL')
LTRE = re.compile('/LT[0-9].')


def spiffify():

    if len(sys.argv) == 1:
        print usage
        sys.exit(2)

    try:
        opts, args = getopt.getopt(sys.argv[1:], "0:1:2:3:4:5:6:7:8:t:s")

    except getopt.GetoptError:
        print "Invalid calling sequence"
        print usage
        sys.exit(2)

    LTDict={}
    thickness="15.000"
    solidify=False

    for opt, value in opts:

        # check for color specifiers
        if opt in ['-'+str(line) for line in range(8)]:

            if value.count(":") != 2:
                print "Invalid color for linetype %s." % opt.lstrip("-")
                sys.exit(2)

            LTDict[opt.lstrip("-")] = value.replace(":"," ")

        if opt == "-t":
            thickness = value

        if opt == "-s":
            solidify=True

    myLineWidth="/mylinewidth " + thickness + " def \n"

    
    for psfile in args:

        foundPL = False
        addedMyPL = False

        # open the file
        try: f = open(psfile,"r")
        except IOError:
            print "ERROR: %s does not exist." % psfile
            sys.exit(2)
        

        # read the lines into an array and close it
        lines=f.readlines()
        f.close()

        # reopen the file for writing
        try: f=open(psfile+".spiffy","w")
        except IOError:
            print "ERROR: %s could not be opened for writing." % psfile
            sys.exit(2)

        # loop over the lines and fix what is needed
        for line in lines:

            outputLine = line

            # add mylinewidth
            if userLineWidthRE.match(line):
                outputLine = line + myLineWidth

            # add myPL
            if not addedMyPL:

                # first, find the PL line
                if not foundPL:
                    if PLRE.match(line):
                        myPLLine = \
                            line.replace("PL","myPL").replace("user","my")

                        foundPL = True

                        outputLine = line

                # we have found and output the PL line.  now we grab the second
                # line in the /PL definition and create the /myPL definition
                else:
                    outputLine = line +  myPLLine + line

                    addedMyPL = True

            # fix LTs
            if LTRE.match(line):
                # use myPL
                outputLine = line.replace("PL","myPL")

                # apply specified colors
                # they seem to always be specified between the 
                # "[...] ... DL" part of the LT definition
                for line, color in LTDict.iteritems():

                    if "/LT"+line in outputLine:
                        irbracket=outputLine.find(']')
                        iD=outputLine.find('D')
                        outputLine = outputLine[:irbracket+1] + " " + \
                            color + " " + outputLine[iD:]

                # make all lines to be a solid linetype
                # that is, remove everything between "[" and "]"
                if solidify:
                    ilbracket=outputLine.find('[')
                    irbracket=outputLine.find(']')
                    outputLine = outputLine[:ilbracket+1] + \
                        outputLine[irbracket:]
                            

            f.write(outputLine)

        f.close()


if __name__ == "__main__":
    spiffify()
