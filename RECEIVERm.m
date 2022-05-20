%RECEIVER MAIN FILE (CONFIGURATION AND PROCESS ROUTINE)
clearvars
close all
addpath(".\..\")

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: SETTINGS ACQUITION        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DAVIDE
inout = InOutInterface();
%Define in-out directories
inout.configCreateSettings("inData");
inout.configSaveResults("outData");
%Read settings
inout.createSettings("in0.json", "PRNpattern.json");
txSymbolRate = inout.settings.chipRate / inout.settings.nChip_x_PRN / inout.settings.nPRN_x_Symbol;

%ToAccess settings : inout.settings, inout.results
% InOutInterface grants that the property ALWAYS contains 
% a valid configuration

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: CONFIGURATION             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%MATTIA
reader = BinaryReader();
%Define input file
reader.configReadFile("binData/testSignals", "T_tracking_1a.bin", inout.settings.quantizationBits);

%GABRIELE
downFilter = DownconverterFilter();
%Define downconverter
downFilter.configDownConverter(inout.settings.fSampling);
%Define filter parameters
% filPassbandStopbandRatio = 1.1;
filRipple_dB = 1;
filAttenuation_dB_dec = 250;
downFilter.configFilter(filRipple_dB,filAttenuation_dB_dec);

%MARCELLO
correlator = CorrelationManager();
correlator.configCorrelatorMatrix(inout.settings.fSampling, inout.settings.nPRN_x_Symbol, ...
                                  inout.settings.nChip_x_PRN, inout.PRNcode, ...
                                  inout.settings.chipRate, inout.settings.SYNCpattern);
%TODO vedere cosa serve per interpolation

%LORENZO
demodulator = Demodulator();
%Define demodulator
demodulator.configCorrelatorValues(inout.settings.fSampling, inout.settings.nPRN_x_Symbol, ...
                                   inout.PRNcode);
%TODO vedere cosa serve per interpolation
%Define packet analysis
demodulator.configMessageAnalyzer(inout.settings.CRCpolynomial,inout.settings.SV_PRN_ID, ...
                                  inout.settings.SYNCpattern, inout.settings.SVIDlength, ...
                                  inout.settings.MIDlength,inout.settings.MBODYlength_ACK, ...
                                  inout.settings.CRClength, inout.settings.TAILpattern);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3: SIGNAL ACQUISITION        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lastSample = min(reader.totSamples, ...
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
smallMaxDoppler = inout.settings.maxDoppler / 2e2; %10
fDopplerResolution = smallMaxDoppler / 40; %1e3
dimCorrMatrix = 1e2; %dimension of output matrices (NOT AFFECTING THE RESULT PRECISION)
thresholdSTD = 3;

while currentSample < lastSample
    %MATTIA
    reader.readFile(currentSample, enlargeForDopplerShift * windowSize);%read file piece and format
    %if currentSample+windowSize>reader.nSamples -> adattare il currentSample
    %->garantire pezzo letto sia lungo windowSize campioni

    % access read data as reader.IQsamples (int_nbit)
    % access as reader.IQsamples_float (single 16bit)
    
    %GABRIELE
    filterBand = txSymbolRate + inout.settings.maxDoppler;
    interFrequency = 0;
    delayLO = 0;
    phase = 0;
    downFilter.downConverter(reader, interFrequency, delayLO, phase);
    downFilter.downFilter(reader, filterBand, inout.settings.chipRate);

    %TEMP BYPASS TILL CODE IS CORRECTED
    correlator.fDoppler = 40;
    correlator.startingSample = 0;
    correlator.initialPhase = 0;
%     break;
    %MARCELLO
    [maxMatrix, meanMatrix, squareMatrix] = ...
                correlator.calcCorrelationMatrix(reader.IQsamples_float, dimCorrMatrix, ...
                                        smallMaxDoppler, fDopplerResolution, currentSample);
    correlator.searchResults
    %access search results as correlator.searchResults
    if correlator.ifAcquiredFindCorrelationPeak(thresholdSTD, inout)
        break;
        %save in inout.results the doopler and the starting time

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

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 4: SIGNAL TRACKING           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lastSymbol = inout.settings.SVIDlength + inout.settings.MIDlength + ...
             inout.settings.MBODYlength_ACK + inout.settings.CRClength + ...
             length(inout.settings.SYNCpattern) + length(inout.settings.TAILpattern);
%segment parameters
currentSymbol = 0;
segmentSize = 1; %number of symbols analyzed togheter
enlargeForDopplerShift = 1.1; %acquire more samples to handle longer symbols
%demodulation parameters
chipFraction = 0.3; %fraction of code shift per tracking
fFraction = 1e-5 / segmentSize; %fraction of doppler frequency shift per tracking
                  %since chiprate is 1Mhz -> fFraction=1e-6 leads to a
                  %doppler shift of ~+-1Hz = 1Mhz*1e-6
coherenceFraction = 0.25; %fraction of symbol period per incoherent detection
all_decodedSymbols = zeros(lastSymbol, 1);

while currentSymbol < lastSymbol
    reader.readFile(correlator.startingSample, ...
                    enlargeForDopplerShift * segmentSize * correlator.nSamples_x_symbolPeriod);%read file piece and format
    
    %GABRIELE
    %here we have segmentSize symbols modulated -> DYNAMIC filter & downcconvertion
    %filterBand = symbolRate + correlator.e_doppler;
    filterBand = 0.6 / correlator.chipPeriod;
    downFilter.downConverter(reader, correlator.fDoppler, ...
                             correlator.startingTime, correlator.initialPhase);
    downFilter.downFilter(reader, filterBand, 1/correlator.chipPeriod);
    %pspectrum(reader.IQsamples_float(:,1),"spectrogram")

    %LORENZO
    %create shifts for delay and doppler
    shifts_delayPRN = int32(correlator.nSamples_x_chipPeriod * chipFraction * [-3 -1 0 1 3]);
    multFactor = 1 + fFraction * [-4 -1 0 1 4]';
    shifts_nSamples_x_symbolPeriod = correlator.nSamples_x_symbolPeriod * multFactor;
    shifts_nSamples_x_chipPeriod = correlator.nSamples_x_chipPeriod * multFactor;
    [decSymbols, idTimeShift, idFreqShift] = ...
                     demodulator.decodeOptimumShift(reader.IQsamples_float, segmentSize, ...
                                                    shifts_delayPRN, shifts_nSamples_x_symbolPeriod, ...
                                                    shifts_nSamples_x_chipPeriod, coherenceFraction);
    decodedSymbols(currentSymbol + 1:currentSymbol + segmentSize) = decSymbols;

    %DAVIDE
    %update the correlation peak with new doppler and delay estimations
    new_samplesChipPeriod = shifts_nSamples_x_chipPeriod(idFreqShift);    
    advancement_startingSample = segmentSize * shifts_nSamples_x_symbolPeriod(idFreqShift) + ...
                                 shifts_delayPRN(idTimeShift);
    correlator.updateCorrelationPeak(new_samplesChipPeriod,advancement_startingSample);
    
    %segment advancement4
    currentSymbol = currentSymbol + segmentSize;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 5: EXTRACTION MESSAGE INFO   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%LORENZO
demodulator.analyzeMessage(decodedSymbols, inout);
%save in inout.results the message splits, the ACKedMessage and ACK flags

%DAVIDE
saved = inout.saveResults("test_results.json");
sprintf("Results saved in %s", saved)

