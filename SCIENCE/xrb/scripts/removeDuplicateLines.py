#!/usr/bin/python
#
"""
  This routine removes the duplicate lines in a text file, inputFile, and
  preserves the order of the lines.

  It turns out that the uniq shell command is capable of doing this
  and as such this routine is just a wrapper for that command.
"""

def removeDuplicateLines(inputFile):
    import sys, os

    # make sure the file is openable
    try: open(inputFile, 'r')
    except IOError:
        print "ERROR: input file %s does not exist." % inputFile
        sys.exit(2)
        
    tempFile = inputFile + '.temp'

    uniqCommand = "uniq -u " + inputFile + " > " + tempFile

    os.system(uniqCommand)

    os.rename(tempFile,inputFile)

if __name__ == "__main__":
    import sys
    
    usage = """
               Calling sequence:
                      %s <inputFile>
            """ % sys.argv[0]

    if len(sys.argv) == 1:
        print usage
        sys.exit()
    else:
        inputFile  = sys.argv[1]
        removeDuplicateLines(inputFile)
