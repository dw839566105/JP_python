# Convert ROOT file to HDF5 file

import numpy as np
import ROOT
import sys
import os
import tables

windowSize = 1029

# Define the database columns
class WaveformData(tables.IsDescription):
    EventID    = tables.Int64Col(pos=0)
    ChannelID  = tables.Int16Col(pos=1)
    Waveform   = tables.Col.from_type('int16', shape=windowSize, pos=2)

class TriggerInfoData(tables.IsDescription):
    EventID    = tables.Int64Col(pos=0)
    Sec        = tables.Int32Col(pos=1)
    NanoSec    = tables.Int32Col(pos=2)

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
TriggerInfoTable = h5file.create_table(group, "TriggerInfo", TriggerInfoData, "Trigger info")
triggerinfo = TriggerInfoTable.row

WaveformTable = h5file.create_table(group, "Waveform", WaveformData, "Waveform")
waveform = WaveformTable.row

GroundTruthTable = h5file.create_table(group, "GroundTruth", GroundTruthData, "GroundTruth")
groundtruth = GroundTruthTable.row

# Loop for ROOT files. 
while True:
    t = ROOT.TChain("Readout")
    tTruth = ROOT.TChain("SimTriggerInfo")
    if FileNo==0:
        filename = baseFileName
    else:
        filename = baseFileName[0:-5]+"_%d.root"%(FileNo)
    if  not os.path.exists(filename):
        break
    print(filename)
    t.Add(filename)
    tTruth.Add(filename)
    print("Processing File %d "%(FileNo)+" ...")
    FileNo = FileNo+1
    
    # Loop for event
    for event in t:
        nChannel = event.ChannelId.size()
        wavSize = event.Waveform.size()
        triggerinfo['EventID'] = event.TriggerNo
        triggerinfo['Sec'] = event.Sec
        triggerinfo['NanoSec'] = event.NanoSec
        triggerinfo.append()
        for i in range(nChannel):
            waveform['EventID'] = event.TriggerNo
            waveform['ChannelID'] = event.ChannelId[i]
            waveform['Waveform'] = event.Waveform[i*windowSize:(i+1)*windowSize]
            waveform.append()
    for event in tTruth:
        for PE in event.PEList:
            groundtruth['EventID'] =  event.TriggerNo
            groundtruth['ChannelID'] =  PE.PMTId
            groundtruth['PETime'] =  PE.HitPosInWindow+10   # Add a 10ns offset. The PE time now is more closer to the waveform. 
            groundtruth.append()

    # Flush into the output file
    TriggerInfoTable.flush()
    WaveformTable.flush()
    GroundTruthTable.flush()

h5file.close()
