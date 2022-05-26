close all

%% SIGNAL GENERATION
clearvars
DEBUG = true;
powerReduce = 6.6
dopplerStart = 45.52
dopplerEnd = dopplerStart + 0.5
signalDirectory = "binData/tempTest"
testSignalName =  "T_tempTest"
TESTPARS = struct("name",testSignalName, ...
                  "dstart",dopplerStart, ...
                  "dend",dopplerEnd, ...
                  "deveryChip",true, ...
                  "dmode",1, ...
                  "envelope",1, ...
                  "envelopePhase",0,...
                  "inNoise",1, ...
                  "outNoise",1, ...
                  "outNoise_Length",1223, ...
                  "outBits_Length",7, ...
                  "outPRN_Length",45, ...
                  "powerReduce",powerReduce);
generateTestSignal

%% TEST RUN
clearvars
DEBUG = true;
outDirectory = "outData/tempTest"
settingsName = "in0.json"
resultsName = "test000000.json"
signalDirectory = "binData/tempTest"
testSignalName =  "T_tempTest.bin"

filterBandMultiplier = 1
reduceMaxDoppler = 100
nTestedDoppler = 100
thresholdSTD = 3
ppSegmentSize = 10
nCoherentFractions = 1

RECEIVER

%clearvars