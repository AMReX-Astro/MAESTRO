#!/usr/bin/python
#
"""
  This routine removes the duplicate lines in a text file, inputFile, and
  preserves the order of the lines.  
"""

def removeDuplicateLines(inputFile, outputFile):
    import sys

    # open the input file if we can
    try: fhIn = open(inputFile, 'r')
    except IOError:
        print "ERROR: input file %s does not exist." % inputFile
        sys.exit(2)

    # create a list to store the non-duplicated lines
    outputLines = []

    # loop over the lines in the file
    for line in fhIn:
        
        # if this line is already in the outputLines list, then we skip it
        # if the file is large, this is going to take a while...
        if line in outputLines: continue

        # otherwise, add it to the list
        outputLines.append(line)

    # close the input file
    fhIn.close()

    # open the output file if we can
    try: fhOut = open(outputFile, 'w')
    except IOError:
        print "ERROR: output file %s could not be created." % outputFile
        sys.exit(2)

    # dump the data
    for line in outputLines:
        fhOut.write(line)

    # close the output file
    fhOut.close()





if __name__ == "__main__":
    import sys
    
    usage = """
               Calling sequence:
                      %s <inputFile> <outputFile>
            """ % sys.argv[0]

    if len(sys.argv) == 1:
        print usage
        sys.exit()
    else:
        inputFile  = sys.argv[1]
        outputFile = sys.argv[2]
        removeDuplicateLines(inputFile,outputFile)
