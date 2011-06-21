#!/usr/bin/env python
import sys

# tex format stuff
header=r"""
\section{Runtime Parameters}

Table~\ref{table:runtime} lists the runtime parameters available to
all MAESTRO problems.  Problem-specific runtime parameters are not
shown here.

%%%%%%%%%%%%%%%%
% symbol table
%%%%%%%%%%%%%%%%

\begin{landscape}

{\small

\renewcommand{\arraystretch}{1.5}
%
\begin{center}
\begin{longtable}{|l|p{5.25in}|l|}
\caption[runtime parameters]{runtime parameters.} \label{table:runtime} \\
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



def make_tex_table():

    # open the file
    try: f = open(paramFile, "r")
    except IOError:
        print "ERROR: %s does not exist" % paramFile
        sys.exit(2)

    # local storage for the parameters
    paramsDict={}
    paramsList=[]
    descr=r""

    # read in the file
    # skip all lines before the first empty line
    foundFirstParam = False

    for line in f:
        if not foundFirstParam:
            if line.isspace():
                # this is the first empty line and we begin reading the file 
                # from here on out
                foundFirstParam = True
                continue
            # no blank line found yet, keep going
            continue

        # land here once we have found the first parameter
        currentParam = Parameter()
        
        # skip blank lines
        if line.isspace(): 
            continue

        elif line.startswith("#"):
            # handle descriptions here
            descr+=line[1:].rstrip().replace("@@",r"\newline")
            continue

        else:

            lineList = line.split()

            currentParam.var=lineList[0]
            currentParam.default=lineList[2].replace("_","\_")
            currentParam.description=descr

            descr=r""

        paramsDict[currentParam.var]=currentParam
        paramsList.append(currentParam.var)

    
    # dump the header
    print header

    # sort the parameters and dump them in latex-fashion
    odd = 1
    for param in sorted(paramsList):
        if (odd == 1):
            print "\\rowcolor{tableShade}"
            odd = 0
        else:
            odd = 1

        print "\\verb= ", \
            paramsDict[param].var, \
            " = & ", \
            paramsDict[param].description, \
            " & ", \
            paramsDict[param].default, \
            r"\\"

    # dump the footer
    print footer
                

if __name__ == "__main__":
    make_tex_table()
