# a simple script to look at some of the parameters and variables that
# yt knows about Maestro data

from yt.mods import *

pf = load("plt23437")


# fields
print "fields:"
for p in  pf.field_list:
    print "  " + `p`

# derived fields
print " "
print "derived fields:"
for p in  pf.derived_field_list:
    print "  " + `p`


# runtime parameters (from the job_info file)
print " "
print "runtime parameters: "
keys = pf.parameters.keys()
for k in sorted(keys):
    print "  {:32}: {}".format(k, pf.parameters[k])

# periodicity?
print " "
print "periodicity: "
print "  ", pf.periodicity




