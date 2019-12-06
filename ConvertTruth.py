# Convert ROOT file to HDF5 file

import numpy as np
import ROOT
import sys
import os
import tables

# Define the database columns

class GroundTruthData(tables.IsDescription):
    EventID    = tables.Int64Col(pos=0)
    ChannelID  = tables.Int16Col(pos=1)
    PETime     = tables.Int16Col(pos=2)

# Automatically add multiple root files created a program with max tree size limitation.
if len(sys.argv)!=3:
    print("Wront arguments!")
    print("Usage: python ConvertTruth.py MCFileName outputFileName")
    sys.exit(1)

baseFileName = sys.argv[1]
outputFileName = sys.argv[2]

ROOT.PyConfig.IgnoreCommandLineOptions = True

FileNo = 0

# Create the output file and the group
h5file = tables.open_file(outputFileName, mode="w", title="OneTonDetector",
                          filters = tables.Filters(complevel=9))
group = "/"

# Create tables

GroundTruthTable = h5file.create_table(group, "GroundTruth", GroundTruthData, "GroundTruth")
groundtruth = GroundTruthTable.row
# Loop for ROOT files. 
count = 0
t = ROOT.TChain("Readout")
tTruth = ROOT.TChain("SimTriggerInfo")
tTruth.Add(baseFileName)
t.Add(baseFileName)
# Loop for event
for event in tTruth:
    count = count + 1
    print(count)
    for PE in event.PEList:
        groundtruth['EventID'] =  event.TriggerNo
        groundtruth['ChannelID'] =  PE.PMTId
        groundtruth['PETime'] =  PE.HitPosInWindow
        groundtruth.append()

# Flush into the output file
GroundTruthTable.flush()

h5file.close()
