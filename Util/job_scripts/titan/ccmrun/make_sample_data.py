import os

prefix = "xrb_3d_w30.72m_plt"


for n in range(20):
    os.mkdir("{}{:05d}".format(prefix, n))

