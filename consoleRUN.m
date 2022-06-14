close all

%% SIGNAL GENERATION
clearvars
DEBUG = true;
powerReduce = 5
dopplerStart = 45.52
dopplerEnd = dopplerStart + 5
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
pause(0.5)
%pause

%% RECEIVER SIMULATION
clearvars
DEBUG = true;
outDirectory = "outData/tempTest"
settingsName = "settings.ini"
resultsName = "test000000.json"
signalDirectory = "binData/tempTest"
testSignalName =  "T_tempTest.bin"

%preciseInterval = 281650
filterBandMultiplier = 1
reducedMaxDoppler = 100 %comment to search on the full interval...
nTestedDoppler = 100 %... but test more frequencies to keep the same resolution
thresholdSTD = 4 %correct threshold for tracking = 3/4
ppSegmentSize = 5
nCoherentFractions = 1

if filterBandMultiplier == 1
    filAttenuation_dB_dec = 390;
elseif filterBandMultiplier == 3
    filAttenuation_dB_dec = 2900;
end

RECEIVER

%clearvars
