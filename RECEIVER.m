%
% First implementation: Mattia Piana & Davide Tomasella
% Review and Testing:
%

%%
%RECEIVER MAIN FILE (CONFIGURATION AND PROCESS ROUTINE)
%DT$ removed clearvars to hande input arguments
%clearvars
close all
addpath(".\..\")

%%
% DEBUG SETTINGS
if ~exist('DEBUG',"var")
    DEBUG = false;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: SETTINGS ACQUITION        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DAVIDE
if ~exist('newResult',"var")
    newResult = true;
end
inout = InOutInterface(DEBUG, newResult);

if ~exist('settingsName',"var")
    settingsName = "in0.json";
end
if ~exist('resultsName',"var")
    resultsName = "test_results.json";
end
if ~exist('outDirectory',"var")
    outDirectory = "outData";
end
%Define in-out directories
inout.configCreateSettings("inData");
inout.configSaveResults(outDirectory);
%Read settings
inout.createSettings(settingsName, "PRNpattern.json");
txSymbolRate = inout.settings.chipRate / inout.settings.nChip_x_PRN / inout.settings.nPRN_x_Symbol;

%ToAccess settings : inout.settings, inout.results
% InOutInterface grants that the property ALWAYS contains 
% a valid configuration

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: CONFIGURATION             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%MATTIA
signal = SignalManager(DEBUG);
if ~exist('testSignalName',"var")
    testSignalName = "T_tracking_1a.bin";
end
if ~exist('signalDirectory',"var")
    signalDirectory = "binData/testSignals";
end
%Define input file
signal.configReadFile(signalDirectory, testSignalName, inout.settings.quantizationBits);

%GABRIELE
downFilter = DownconverterFilter(DEBUG);
%Define downconverter
downFilter.configDownConverter(inout.settings.fSampling);
%Define filter parameters
% filPassbandStopbandRatio = 1.1;
filRipple_dB = 1;
filAttenuation_dB_dec = 250;
if ~exist('filterBandMultiplier',"var")
    filterBandMultiplier = 1;
end
downFilter.configFilter(filRipple_dB, filAttenuation_dB_dec);

%MARCELLO
correlator = CorrelationManager(DEBUG);
correlator.configCorrelatorMatrix(inout.settings.fSampling, inout.settings.nPRN_x_Symbol, ...
                                  inout.settings.nChip_x_PRN, inout.PRNcode, ...
                                  inout.settings.chipRate, inout.settings.SYNCpattern, ...
                                  inout.settings.quantizationBits);

%LORENZO
%Define demodulator
tracker = TrackingManager(DEBUG);
MMSEalpha = 0.9;
tracker.configCorrelatorValues(inout.settings.fSampling, inout.settings.nPRN_x_Symbol, ...
                               inout.PRNcode, MMSEalpha, inout.settings.quantizationBits);
%Define packet analysis
analyzer = MessageAnalyzer(DEBUG);
analyzer.configMessageAnalyzer(inout.settings.CRCpolynomial, inout.settings.SV_PRN_ID, ...
                               inout.settings.SYNCpattern, inout.settings.SVIDlength, ...
                               inout.settings.MIDlength, inout.settings.MBODYlength_ACK, ...
                               inout.settings.CRClength, inout.settings.TAILpattern);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3: SIGNAL ACQUISITION        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lastSample = min(signal.totSamples, ...
                 inout.settings.scenarioDuration * inout.settings.fSampling);
%window parameters
currentSample = 0;
nSamples_x_symbolPeriod = (1 / txSymbolRate) * inout.settings.fSampling;
nSymbols_per_sync = length(inout.settings.SYNCpattern);
kLength = 2;
windowSize = kLength * nSamples_x_symbolPeriod * nSymbols_per_sync;
enlargeForDopplerShift = 1.1; %acquire more samples to handle longer symbols
windowsAdvancement = windowSize - nSamples_x_symbolPeriod * nSymbols_per_sync; %+1

%acquisition paramters
if ~exist('thresholdSTD',"var")
    thresholdSTD = 3;
end
if ~exist('reducedMaxDoppler',"var")
    reducedMaxDoppler = inout.settings.maxDoppler / 1e2; %(100kHz auto-reduced)
end
if ~exist('centralDoppler',"var")
    centralDoppler = 0;
end
if ~exist('nTestedDoppler',"var")
    nTestedDoppler = 100;
end
fDopplerResolution = 2* reducedMaxDoppler / nTestedDoppler; %1e3
dimCorrMatrix = 100; %dimension of output matrices (NOT AFFECTING THE RESULT PRECISION)

while currentSample < lastSample
    %MATTIA
    signal.readFile(currentSample, enlargeForDopplerShift * windowSize);%read file piece and format
    % access read data as reader.IQsamples (int_nbit)
    % access as reader.IQsamples_float (single 16bit)
    
    %GABRIELE
    %No downconversion because the signals are already in base-band
    interFrequency = 0;
    delayLO = 0;
    phase = 0;
    filterBand = filterBandMultiplier * (txSymbolRate + inout.settings.maxDoppler);
    %downFilter.downConverter(signal, interFrequency, delayLO, phase);
    downFilter.downFilter(signal, filterBand, inout.settings.chipRate);

    %GABRIELE
    [maxMatrix, meanMatrix, squareMatrix] = ...
                correlator.calcCorrelationMatrix(signal.IQsamples_float, dimCorrMatrix, ...
                                                 reducedMaxDoppler, centralDoppler, ...
                                                 fDopplerResolution, currentSample);
    %access search results as correlator.searchResults
    correlator.searchResults

    %DAVIDE
    if correlator.ifAcquiredFindCorrelationPeak(thresholdSTD, inout)
        %save in inout.results the doopler and the starting time
        break;

        % access optimal parameters as correlator.fDoppler,
        %                              correlator.symbolPeriod,
        %                              correlator.chipPeriod,
        %                              correlator.nSamples_x_symbolPeriod
        %                              correlator.nSamples_x_chipPeriod
        %                              correlator.startingSample,
        %                              correlator.startingTime
    end
    
    %window advancement
    currentSample = currentSample + windowsAdvancement;
end

if  isempty(correlator.startingSample) || ...
            correlator.startingSample > lastSample - windowSize
    inout.results.ACQUISITION_OK = false;
    warning("NO SIGNALS FOUND")
    return
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 4: SIGNAL TRACKING           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lastSymbol = inout.settings.SVIDlength + inout.settings.MIDlength + ...
             inout.settings.MBODYlength_ACK + inout.settings.CRClength + ...
             length(inout.settings.SYNCpattern) + length(inout.settings.TAILpattern);
%segment parameters
if ~exist('ppSegmentSize',"var")
    ppSegmentSize = 10; %number of symbols analyzed togheter
end
currentSymbol = 0;
enlargeForDopplerShift = 1.1; %acquire more samples to handle longer symbols

%demodulation parameters
if ~exist('nCoherentFractions',"var")
    nCoherentFractions = 1; %fraction of symbol period per incoherent detection
end
chipFraction = 0.1; %fraction of code shift per tracking
fFraction = 1e-6 / sqrt(ppSegmentSize); %fraction of doppler frequency shift per tracking
                                      %since chiprate is 1Mhz -> fFraction=1e-6 leads to a
                                      %doppler shift of ~+-1Hz = 1Mhz*1e-6
plotVector1 = [-20 -8 -4 -2 0 2 4 8 20];
plotVector2 = [-40 -20 -10 -5 -1 0 1 5 10 20 40]';
delayVector = [-8 -2 0 2 8];
dopplerVector = [-12 -2 0 2 12]';
all_decodedSymbols = zeros(lastSymbol, 1);
inout.results.TRACKING_OK = true; %initialize tracking OK

while currentSymbol < lastSymbol
    signal.readFile(correlator.startingSample, ...
                    enlargeForDopplerShift * ppSegmentSize * correlator.nSamples_x_symbolPeriod);%read file piece and format
    
    %GABRIELE
    %here we have segmentSize symbols modulated -> DYNAMIC filter & downcconvertion
    %filterBand = symbolRate + correlator.e_doppler;
    %DT$ startingTime=0 because we estimate the phase delay thanks to doppler
    filterBand = filterBandMultiplier * correlator.chipPeriod;
    downFilter.downConverter(signal, correlator.fDoppler, ...
                             0, correlator.initialPhase);
    downFilter.downFilter(signal, filterBand, 1/correlator.chipPeriod);
    %plot(signal.IQsamples(:,1))

    %LORENZO
    %create shifts for delay and doppler
    shifts_delayPRN = int32(correlator.nSamples_x_chipPeriod * chipFraction * delayVector);
    multFactor = 1 + fFraction * dopplerVector;
    shifts_nSamples_x_symbolPeriod = correlator.nSamples_x_symbolPeriod * multFactor;
    shifts_nSamples_x_chipPeriod = correlator.nSamples_x_chipPeriod * multFactor;
    %calculate best shift and demodulate
    [decSymbols, idTimeShift, idFreqShift] = ...
                     tracker.decodeOptimumShift(signal.IQsamples_float, ppSegmentSize, ...
                                                shifts_delayPRN, shifts_nSamples_x_symbolPeriod, ...
                                                shifts_nSamples_x_chipPeriod, nCoherentFractions, ...
                                                inout);
    %save demodulated symbols
    decodedSymbols(currentSymbol + 1:currentSymbol + ppSegmentSize) = decSymbols;

    %DAVIDE
    %update the correlation peak with new doppler and delay estimations
    new_samplesChipPeriod = shifts_nSamples_x_chipPeriod(idFreqShift);    
    advancement_startingSample = ppSegmentSize * shifts_nSamples_x_symbolPeriod(idFreqShift) + ...
                                 shifts_delayPRN(idTimeShift);
    correlator.updateCorrelationPeak(new_samplesChipPeriod, ...
                                     advancement_startingSample);
    
    %segment advancement
    currentSymbol = currentSymbol + ppSegmentSize;
end

inout.results.estimatedDopplerEnd = round(correlator.fDoppler, 6, "significant");
%save in inout.results the doppler at the end of the message

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 5: EXTRACTION MESSAGE INFO   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%LORENZO
demodOK = analyzer.analyzeMessage(decodedSymbols, inout);
sprintf("DEMODULATION RESULT %d", demodOK)

%save in inout.results the message splits, the ACKedMessage and ACK flags
inout.results

%DAVIDE
saved = inout.saveResults(resultsName);
sprintf("Results saved in %s", saved)

