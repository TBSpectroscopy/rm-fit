#---------------------------------------------------------------------------------------------------
# Thibault BERTIN
# Atomic and Molecular Physics
# Center for Astrophysics | Harvard & Smithsonian
# 60 Garden St., 02138 MA, USA
# E-mail: thibault.bertin@cfa.harvard.edu
#---------------------------------------------------------------------------------------------------

import os
import sys
import datetime
import src.msfp as msfp

print("RM-Fit")
if len(sys.argv) == 1 or sys.argv[1] in ["-h", "--help"]:
    print("Usage: python rm-fit.py [OPTION] [FILE]\n"\
            "Options:\n"\
            "   -h, --help      Print this help text and exit\n"\
            "   -c              Calculate spectra using given parameters\n"\
            "   -f              Fit theoritical curves to observed spectra")
    sys.exit()

print("Run time: {}".format(datetime.datetime.now().isoformat(sep = " ", timespec="minutes")))
if len(sys.argv) == 2 :
    print("ERROR - Missing argument: option and  path to an input file must be provided\n\nTry \"python rm-fit -h\" or \"python rm-fit --help\" for more information")
    sys.exit()
elif sys.argv[1] not in ["-c", "-f"]:
    print("ERROR: unrecognized option \"{}\"\nTry \"python rm-fit -h\" or \"python rm-fit --help\" for more information".format(sys.argv[1]))
    sys.exit()

inpdir = ".{}{}".format(os.sep, os.path.dirname(sys.argv[2]))
inppath = os.path.basename(sys.argv[2])
os.chdir(inpdir)
msfp.msfp(inppath, sys.argv[1])
