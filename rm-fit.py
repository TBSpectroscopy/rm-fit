import os
import sys
import datetime
import src.msfp as msfp

print("msfp.py")
if len(sys.argv) == 1 or sys.argv[1] == "-h" or sys.argv[1] == "--help":
    print("Usage: msfp [OPTION] [FILE]\n"\
            "Options:\n"\
            "   -h, --help      Print this help text and exit\n"\
            "   -c              Calculate spectra using given parameters\n"\
            "   -f              Fit theoritical curves to observed spectra")
    sys.exit()

print("Run time: {}".format(datetime.datetime.now().isoformat(sep = " ", timespec="minutes")))
if len(sys.argv) == 2 :
    print("ERROR - Missing argument: the path (relative to the directory containing " \
        "\"msfp.py\") to a \"control file\" must be provided")
    sys.exit()
inpdir = os.path.dirname(sys.argv[2])
inppath = os.path.basename(sys.argv[2])
os.chdir(inpdir)
msfp.msfp(inppath, sys.argv[1])
