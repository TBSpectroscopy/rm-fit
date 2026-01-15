#-------------------------------------------------------------------------------------------------
#
# Thibault BERTIN
# Atomic and Molecular Physics
# Center for Astrophysics | Harvard & Smithsonian
# 60 Garden St., 02138 MA, USA
# E-mail: thibault.bertin@cfa.harvard.edu
#
#-------------------------------------------------------------------------------------------------


import sys
import os
from inspect import getsourcefile


src_dir = os.path.dirname(os.path.abspath(getsourcefile(lambda: 0)))
sys.path.append(src_dir)
