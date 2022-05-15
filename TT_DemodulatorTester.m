%% configs base
clearvars
close all
addpath(".\..\")

inout = InOutInterface();
%Define in-out directories
inout.configCreateSettings("inData");
inout.configSaveResults("outData");
%Read settings
inout.createSettings("in0.json", "PRNpattern.json");
txSymbolRate = inout.settings.chipRate / inout.settings.nChip_x_PRN / inout.settings.nPRN_x_Symbol;

reader = BinaryReader();
%Define input file
reader.configReadFile("binData/testSignals", "T_tracking_1_220515_114108.bin", inout.settings.quantizationBits);

%% configs actions
downFilter = DownconverterFilter();
%Define downconverter
downFilter.configDownConverter(inout.settings.fSampling);
%Define filter parameters
filPassbandStopbandRatio = 1.1;
filRipple_dB = 0.3;
filAttenuation_dB = 30;
downFilter.configFilter(filPassbandStopbandRatio, filRipple_dB, filAttenuation_dB);

correlator = Correlator();
correlator.configCorrelatorMatrix(inout.settings.fSampling, inout.settings.nPRN_x_Symbol, ...
                                  inout.settings.nChip_x_PRN, inout.PRNcode, ...
                                  inout.settings.chipRate, inout.settings.SYNCpattern);

demodulator = TT_Demodulator();
%Define demodulator
demodulator.configCorrelatorValues(inout.settings.fSampling, inout.settings.nPRN_x_Symbol, ...
                                   inout.PRNcode);
%TODO vedere cosa serve per interpolation
%Define packet analysis
demodulator.configMessageAnalyzer(inout.settings.CRCpolynomial,inout.settings.SV_PRN_ID, ...
                                  inout.settings.SYNCpattern, inout.settings.SVIDlength, ...
                                  inout.settings.MIDlength,inout.settings.MBODYlength_ACK, ...
                                  inout.settings.CRClength, inout.settings.TAILpattern);

%% setup parameter
%correlator bypass
correlator.fDoppler = 15.23;
correlator.startingSample = 0;

%T_tracking_1/6 contains 20 symbols
lastSymbol = 20;

%segment parameters
currentSymbol = 0;
segmentSize = 3; %number of symbols analyzed togheter
enlargeForDopplerShift = 1.1; %acquire more samples to handle longer symbols
%demodulation parameters
chipFraction = 0.5; %fraction of code shift per tracking
fFraction = 1e-6; %fraction of doppler frequency shift per tracking
                  %since chiprate is 1Mhz -> fFraction=1e-6 leads to a
                  %doppler shift of ~+-1Hz = 1Mhz*1e-6
coherenceFraction = 1; %fraction of symbol period per incoherent detection
all_decodedSymbols = zeros(lastSymbol, 1);

%% read and pre-process signal
%read file
reader.readFile(correlator.startingSample, ...
                    enlargeForDopplerShift * segmentSize * correlator.nSamples_x_symbolPeriod);%read file piece and format

filterBand = 1.01 / correlator.chipPeriod;
%no downconvertion for T_tracking_1/4
if false
    downFilter.downConverter(reader, correlator.fDoppler, ...
                                               correlator.startingTime);
    downFilter.downFilter(reader, filterBand, 1/correlator.chipPeriod);
end
plot(reader.IQsamples(:,1))

%% demodulate
shifts_delayPRN = int32(correlator.nSamples_x_chipPeriod * [-2 * chipFraction , -chipFraction, ...
                                                            0, chipFraction, 2 * chipFraction]);
%use constant frequency
if true
    shifts_nSamples_x_symbolPeriod = correlator.nSamples_x_symbolPeriod * [1 - fFraction, 1, 1 + fFraction]';
    shifts_nSamples_x_chipPeriod = correlator.nSamples_x_chipPeriod * [1 - fFraction, 1, 1 + fFraction]';
else
    shifts_nSamples_x_symbolPeriod = correlator.nSamples_x_symbolPeriod;
    shifts_nSamples_x_chipPeriod = correlator.nSamples_x_chipPeriod;
end

%% tracking

[all_decSymbols, all_idTimeShift, all_idFreqShift] = ...
                 demodulator.decodeOptimumShift(reader.IQsamples_float, segmentSize, ...
                                                shifts_delayPRN, shifts_nSamples_x_symbolPeriod, ...
                                                shifts_nSamples_x_chipPeriod, coherenceFraction)

%% tracking splitted
inSamples = reader.IQsamples_float;

%create PRN with different length due to different dopplers
PRNsampled = demodulator.getPRNFromChipPeriods(shifts_nSamples_x_chipPeriod, segmentSize);

%adapt the length of the samples to the PRN, horizontal vector
mySamples = demodulator.adaptSamplesToPRNlength(inSamples, size(PRNsampled, 2));

%add different delays to the PRN, sampled PRN with shifts in time and frequency
shiftedPRNsampled = demodulator.createShiftedPRN(PRNsampled, shifts_delayPRN);
%plot(shiftedPRNsampled(1:2,end-1e4:1:end)')

%in-phase & quadrature multicorrelation
corrI = demodulator.normMultiply(shiftedPRNsampled, mySamples(1, :));
corrQ = demodulator.normMultiply(shiftedPRNsampled, mySamples(2, :));

%coherent integration, over the rows there are the coherent sums
coherentCorrI = demodulator.sumOverCoherentFraction(corrI, shifts_nSamples_x_symbolPeriod, ...
                                            coherenceFraction, segmentSize);
coherentCorrQ = demodulator.sumOverCoherentFraction(corrQ, shifts_nSamples_x_symbolPeriod, ...
                                            coherenceFraction, segmentSize);

%non-coherent integration, column vector
noncoherentCorr = sum(coherentCorrI .^ 2 + coherentCorrQ .^ 2, 2);

%find max
[~, idMax] = max(noncoherentCorr,[], 1);           
[idDoppler, idShift] = ind2sub([size(shifts_nSamples_x_chipPeriod,1) size(shifts_delayPRN,2)], idMax);         
%NOTE: DOT product: sommatoria[(Ie-Il)*Ip] - sommatoria[(Qe-Ql)*qp],
%if <0 detector is late, if >0 too early

%complete correlation over symbols
%TODO channel inversion? linear estimator?
bestCorrI = demodulator.sumFractionsOverSymbols(coherentCorrI(idMax, :), coherenceFraction);
bestCorrQ = demodulator.sumFractionsOverSymbols(coherentCorrQ(idMax, :), coherenceFraction);

%decoding
decSymbols = (2 * (bestCorrI > 0) - 1)'; % column vector of decoded symbols +1,-1            


%% update
if idShift~=all_idTimeShift 
    warning('DIFFERENT %s %s',idShift, all_idTimeShift);
    disp([shifts_delayPRN(idShift) shifts_delayPRN(all_idTimeShift)]);
end
if idDoppler~=all_idFreqShift 
    warning('DIFFERENT %s %s',idDoppler, all_idFreqShift);
    disp([shifts_nSamples_x_symbolPeriod(idDoppler) shifts_nSamples_x_symbolPeriod(all_idFreqShift)]);
end
if decSymbols(1)~=all_decSymbols(1) 
    warning('DIFFERENT %s %s',decSymbols, all_decSymbols);
end

all_decodedSymbols(currentSymbol + 1:currentSymbol + segmentSize) = all_decSymbols;

%update the correlation peak with new doppler and delay estimations
new_samplesSymbolPeriod = shifts_nSamples_x_symbolPeriod(all_idFreqShift);    
advancement_startingSample = segmentSize * shifts_nSamples_x_symbolPeriod(all_idFreqShift) + ...
                             shifts_delayPRN(all_idTimeShift);
correlator.updatePeak(new_samplesSymbolPeriod,advancement_startingSample);

