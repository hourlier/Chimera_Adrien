import sys

from larlite import larlite as fmwk

# Create ana_processor instance
my_proc = fmwk.ana_processor()

# Set input root file
# set default file if none is specified
if len(sys.argv) == 1:
    #my_proc.add_input_file("/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/DeepLearning/DLwork/DataFiles/larlite_extBNB9131runs_cosmic_trained_only_on_mc_pscore_0.99_1598evts_23aug2016.root")
    #my_proc.add_input_file("/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/DeepLearning/DLwork/DataFiles/larlite_pandoraNu_20170210_000045_475610.root")
    my_proc.add_input_file("/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/DeepLearning/DLwork/DataFiles/larlite_pandoraNu_20170210_000730_864687.root")

# use specified files
if len(sys.argv) > 1:
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
my_proc.run(0);

# done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)
