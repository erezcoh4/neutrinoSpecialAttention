import sys
from ROOT import gSystem
gSystem.Load("libneutrinoSpecialAttention_nuSpecAttPackage")
from ROOT import sample

try:

    print "PyROOT recognized your class %s" % str(sample)

except NameError:

    print "Failed importing nuSpecAttPackage..."

sys.exit(0)

