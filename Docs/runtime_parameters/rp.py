#!/usr/bin/env python
import os
import sys

# tex format stuff
Mheader=r"""
\label{ch:parameters}

\begin{landscape}
"""

header=r"""
{\small

\renewcommand{\arraystretch}{1.5}
%
\begin{center}
\begin{longtable}{|l|p{5.25in}|l|}
\caption[@@catname@@]{@@catname@@} \label{table: @@catname@@ runtime} \\
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

maestro_top = "../../"

# we list the files to parse here in a tuple, consisting of the display name 
# (for the LaTeX table) and the path
param_paths = [("main Maestro parameters", "Source/"),
               (r"gamma law general  EOS", "Microphysics/EOS/gamma_law_general/"),
               ("generic network parameters", "Microphysics/networks/"),
               (r"constant conductivity", "Microphysics/conductivity/constant")]


class Parameter(object):
    # container class for the parameters
    
    def __init__(self):
        self.var = ""
        self.default = ""
        self.description = []
        self.category = ""

    def value(self):
        """ the value is what we sort based on """
        return self.category + "." + self.var

    def __cmp__(self, other):
        return cmp(self.value(), other.value())

    def __str__(self):
        return "{}".format(self.var)


def make_tex_table(param_info):

    display_name, param_dir = param_info

    param_file = os.path.join(maestro_top, param_dir, "_parameters")

    # open the file
    try: f = open(param_file, "r")
    except IOError:
        print "ERROR: %s does not exist" % param_file
        sys.exit(2)

    # local storage for the parameters
    params_list=[]
    descr=r""
    category=""

    # read in the file
    # skip all lines before the first empty line
    found_first_param = False

    line = f.readline()
    while line:

        # we assume that parameter files have a descriptive heading
        # before the parameters start in those cases, we want to skip
        # everything before the first blank line.  the next comment
        # will then be interpreted either as the category or as the
        # description for the parameter

        if not found_first_param:
            if line.isspace():
                # this is the first empty line and we begin reading the file 
                # from here on out
                found_first_param = True
                line = f.readline()
                continue

            # no blank line found yet, keep going
            line = f.readline()
            continue

        # land here once we have found the first parameter
        current_param = Parameter()
        
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
            if not line.startswith("#------"):
                sys.exit("ERROR: category block not formatted correctly")

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

            current_param.var = lineList[0]
            current_param.default = lineList[2].replace("_","\_")
            current_param.description = descr
            current_param.category = category

            descr=r""

        
        # store the current parameter in the list
        params_list.append(current_param)
        
        line = f.readline()


    # sort the parameters and dump them in latex-fashion.  Group things by category
    current_category = -1
    start = 1

    for param in sorted(params_list):

        if not param.category == current_category:
            if not start == 1:
                print footer

            current_category = param.category
            odd = 1
            if param.category == "":
                cat_header = header.replace("@@catname@@", display_name + " parameters.")
            else:
                cat_header = header.replace("@@catname@@", param.category + " parameters.")
            print cat_header
            start = 0

        if odd == 1:
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

if __name__ == "__main__":

    # dump the main header
    print Mheader

    for p in param_paths:
        make_tex_table(p)

    print Mfooter
