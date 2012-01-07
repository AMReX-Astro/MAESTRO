#!/usr/bin/env python
import sys

# tex format stuff
Mheader=r"""
\section{Runtime Parameters}

Table~\ref{table:runtime} lists the runtime parameters available to
all MAESTRO problems.  Problem-specific runtime parameters are not
shown here.

%%%%%%%%%%%%%%%%
% symbol table
%%%%%%%%%%%%%%%%

\begin{landscape}
"""

header=r"""
{\small

\renewcommand{\arraystretch}{1.5}
%
\begin{center}
\begin{longtable}{|l|p{5.25in}|l|}
\caption[runtime parameters]{@@catname@@} \label{table:runtime} \\
%
\hline \multicolumn{1}{|c|}{\textbf{parameter}} & 
       \multicolumn{1}{ c|}{\textbf{description}} & 
       \multicolumn{1}{ c|}{\textbf{default value}} \\ \hline 
\endfirsthead

\multicolumn{3}{c}%
{{\tablename\ \thetable{}---continued}} \\
\hline \multicolumn{1}{|c|}{\textbf{parameter}} & 
       \multicolumn{1}{ c|}{\textbf{description}} & 
       \multicolumn{1}{ c|}{\textbf{default value}} \\ \hline 
\endhead

\multicolumn{3}{|r|}{{\em continued on next page}} \\ \hline
\endfoot

\hline 
\endlastfoot

"""

footer=r"""

\end{longtable}
\end{center}

} % ends \small
"""

Mfooter=r"""
\end{landscape}

%

"""

fParallel="../../.."
paramFile=fParallel + "/MAESTRO/_parameters"


# container class for the parameters
class Parameter:
    def __init__(self):
        self.var=""
        self.default=""
        self.description=[]
        self.category=""

    def value(self):
        """ the value is what we sort based on """
        return self.category + "." + self.var

    def __cmp__(self, other):
        return cmp(self.value(), other.value())


def make_tex_table():

    # open the file
    try: f = open(paramFile, "r")
    except IOError:
        print "ERROR: %s does not exist" % paramFile
        sys.exit(2)

    # local storage for the parameters
    paramsList=[]
    descr=r""
    category=""

    # read in the file
    # skip all lines before the first empty line
    foundFirstParam = False

    line = f.readline()
    while (line):

        if not foundFirstParam:
            if line.isspace():
                # this is the first empty line and we begin reading the file 
                # from here on out
                foundFirstParam = True
                line = f.readline()
                continue

            # no blank line found yet, keep going
            line = f.readline()
            continue

        # land here once we have found the first parameter
        currentParam = Parameter()
        
        # skip blank lines
        if line.isspace(): 
            line = f.readline()
            continue

    
        # look for category definition
        elif line.startswith("#------"):

            # the next line should be the category definition
            line = f.readline()
            index = line.find(":")
            category = line[index+1:]

            # following this is another #---------
            line = f.readline()
            if (not line.startswith("#------")):
                print "ERROR: category block not formatted correctly"
                sys.exit(2)

            line = f.readline()
            continue

        # find the description
        elif line.startswith("#"):

            # handle descriptions here
            descr+=line[1:].rstrip().replace("@@",r"\newline")
            line = f.readline()
            continue

        else:

            lineList = line.split()

            currentParam.var=lineList[0]
            currentParam.default=lineList[2].replace("_","\_")
            currentParam.description=descr
            currentParam.category=category

            descr=r""

        
        # store the current parameter in the list
        paramsList.append(currentParam)
        
        
        # get the next line
        line = f.readline()

    
    # dump the main header
    print Mheader

    # sort the parameters and dump them in latex-fashion.  Group things by category
    currentCategory = ""
    start = 1

    for param in sorted(paramsList):

        if (not param.category == currentCategory):
            if (not start == 1):
                print footer

            currentCategory = param.category
            odd = 1
            catHeader = header.replace("@@catname@@", param.category + " parameters.")
            print catHeader
            start = 0

        if (odd == 1):
            print "\\rowcolor{tableShade}"
            odd = 0
        else:
            odd = 1

        print "\\verb= ", \
            param.var, \
            " = & ", \
            param.description, \
            " & ", \
            param.default, \
            r"\\"

    # dump the footer
    print footer
    print Mfooter

if __name__ == "__main__":
    make_tex_table()
