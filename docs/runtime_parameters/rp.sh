#!/bin/ksh

# dump the header
cat <<EOF

\section{Runtime Parameters}

%%%%%%%%%%%%%%%%
% symbol table
%%%%%%%%%%%%%%%%

\renewcommand{\arraystretch}{1.5}
%
\begin{center}
\begin{longtable}{|l|p{3.25in}|l|}
\caption[runtime parameters]{runtime parameters.} \label{table:runtime} \\\\
%
\hline \multicolumn{1}{|c|}{\textbf{parameter}} & 
       \multicolumn{1}{ c|}{\textbf{description}} & 
       \multicolumn{1}{ c|}{\textbf{default value}} \\\\ \hline 
\endfirsthead

\multicolumn{3}{c}%
{{\tablename\ \thetable{}---continued}} \\\\
\hline \multicolumn{1}{|c|}{\textbf{parameter}} & 
       \multicolumn{1}{ c|}{\textbf{description}} & 
       \multicolumn{1}{ c|}{\textbf{default value}} \\\\ \hline 
\endhead

\multicolumn{3}{|r|}{{\em continued on next page}} \\\\ \hline
\endfoot

\hline 
\endlastfoot

EOF

# make the table entries
parameters=$(grep namelist probin.f90 | awk '{print $3}' | sort)

for i in ${parameters}
do
    value=$(grep -m 1 "^[ ]*$i =" probin.f90 | cut -d'=' -f2 | cut -d"!" -f1)

    echo "\\verb= " $i "=" " &   & " `echo ${value} | sed 's/\_/\\\_/g'` "\\\\"
done

# dump the end
cat <<EOF

\end{longtable}
\end{center}
%

EOF
