import sys
from ROOT import gSystem
gSystem.Load("libChimera_Adrien_ChimeraProtonTracks")
from ROOT import sample

try:

    print "PyROOT recognized your class %s" % str(sample)

except NameError:

    print "Failed importing ChimeraProtonTracks..."

sys.exit(0)

