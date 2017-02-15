import sys

from larlite import larlite as fmwk

# Create ana_processor instance
my_proc = fmwk.ana_processor()

# Set input root file
# set default file if none is specified
if len(sys.argv) <= 2:
    my_proc.add_input_file("/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/DeepLearning/DLwork/DataFiles/larlite_pandoraNu_20170211_205524_331124.root")

Nevt = -1
if len(sys.argv) == 2:
    Nevt = int(sys.argv[1]);

# use specified files
if len(sys.argv) > 2:
    for x in xrange(len(sys.argv)-1):
        my_proc.add_input_file(sys.argv[x+1])

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kREAD)

# Specify output root file name
my_proc.set_ana_output_file("RESULTS.root");

# Attach an analysis unit ... here we use a base class which does nothing.
# Replace with your analysis unit if you wish.
my_proc.add_process(fmwk.FindProtonTracks())
#my_proc.add_process(fmwk.ana_base())

print
print  "Finished configuring ana_processor. Start event loop!"
print


# Let's run it.
if Nevt <= -1 :
    my_proc.run(0);
if Nevt > 0 :
    my_proc.run(0,Nevt);

# done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)
