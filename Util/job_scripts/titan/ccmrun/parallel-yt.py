#!/usr/bin/env python

import os

num_jobs = 6


prefix = "xrb_3d_w30.72m_plt"


# the command to execute on the file -- leave a {} to indicate
# where the filename gets inserted
vis_command = "python ./vol-perspective-rotate.py {} vz"


# wild card find all the plotfile directories
pfiles = []
pngs = []
for d in os.listdir(os.getcwd()):
    if os.path.isdir(d) and d.startswith(prefix):
        pfiles.append(d)

    if os.path.isfile(d) and d.endswith(".png"):
        pngs.append(d)    


# remove an for which we already have a .png file
pfiles_to_process = []
for d in pfiles:
    found = False
    for image in pngs:
        if d in image:
            found = True
            break
        
    if not found: pfiles_to_process.append(d)
    

# main loop
while len(pfiles_to_process) > 0:

    current = []
    
    # pop off num_job files
    for i in range(min(num_jobs,len(pfiles_to_process))):
        current.append(pfiles_to_process.pop())
        
    # write a script "runtask.sh" that has the command needed to
    # run each of these, ending with "wait"
    print "processing ", current
    
    try: f = open("runtask.sh", "w")
    except:
        sys.exit("ERROR: cannot write the runtask.sh")

    f.write(r"#!/bin/bash" + "\n")
    for c in current:
        f.write(vis_command.format(c) + " &\n")
    f.write("wait\n")
    f.close()
    
    # run with ccmrun
    os.system("ccmrun runtask.sh")

    

    

