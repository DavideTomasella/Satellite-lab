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
reader.configReadFile("binData/testSignals", "T_tracking_1.bin", inout.settings.quantizationBits);
packet = int16([0  1  0  1  1  0  0  0  ...
                0  0  0  0  0  0  0  1  ...
                0  0  1  1  0  1  0  1  ...
                0  0  1  1  0  0  0  1  ...
                1  1  0  0  1  1  0  0  ...
                1  0  1  1  1  0  1  1  ...
                1  0  0  0  0  1  1  0  ...
                1  1  1  1  0  0  1  0  ...
                0  0  1  0  1  0  0  0  ...
                0  1  0  0  0  0  0  0]);

%% configs actions
downFilter = DownconverterFilter();
%Define downconverter
downFilter.configDownConverter(inout.settings.fSampling);
%Define filter parameters
%filPassbandStopbandRatio = 1.1;
filRipple_dB = 1;
filAttenuation_dB = 250;
downFilter.configFilter(filRipple_dB, filAttenuation_dB);

correlator = CorrelationManager();
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
dopplerError = 10;
correlator.fDoppler = 15.23 + dopplerError;
correlator.startingSample = 0;

%T_tracking_1/6 contains 20 symbols
lastSymbol = 80;

%segment parameters
currentSymbol = 0;
segmentSize = 1; %number of symbols analyzed togheter
enlargeForDopplerShift = 1.1; %acquire more samples to handle longer symbols
%demodulation parameters
chipFraction = 0.1; %fraction of code shift per tracking
fFraction = 1e-5 / sqrt(segmentSize); %fraction of doppler frequency shift per tracking
                  %since chiprate is 1Mhz -> fFraction=1e-6 leads to a
                  %doppler shift of ~+-1Hz = 1Mhz*1e-6
coherenceFraction = 1; %fraction of symbol period per incoherent detection
all_decodedSymbols = zeros(lastSymbol, 1);

i=1;
DopplerCorrectionEvolution=zeros(1,lastSymbol/segmentSize);
DelayShiftEvolution=zeros(1,lastSymbol/segmentSize);
FreqShiftEvolution=zeros(1,lastSymbol/segmentSize);

%%
while currentSymbol < lastSymbol
    
DopplerCorrectionEvolution(i)=correlator.fDoppler;
%% read and pre-process signal
%read file
reader.readFile(correlator.startingSample, ...
                    enlargeForDopplerShift * segmentSize * correlator.nSamples_x_symbolPeriod);%read file piece and format

filterBand = 0.6 / correlator.chipPeriod;
%no downconvertion for T_tracking_1/4
if true
    downFilter.downConverter(reader, correlator.fDoppler, ...
                                               correlator.startingTime);
    downFilter.downFilter(reader, filterBand, 1/correlator.chipPeriod);
end
%transform the input in unitary vectors
%reader.IQsamples_float(:,1) = reader.IQsamples_float(:,1) / max(reader.IQsamples_float(:,1));
%plot(reader.IQsamples(:,1))

%% demodulate
shifts_delayPRN = int32(correlator.nSamples_x_chipPeriod * chipFraction * ...
    [-8 -2 0 2 8]);

%use constant frequency
if true
    multFactor = 1 + fFraction * [-5 -1 0 1 5]';
    shifts_nSamples_x_symbolPeriod = correlator.nSamples_x_symbolPeriod * multFactor;
    shifts_nSamples_x_chipPeriod = correlator.nSamples_x_chipPeriod * multFactor;
else
    shifts_nSamples_x_symbolPeriod = correlator.nSamples_x_symbolPeriod;
    shifts_nSamples_x_chipPeriod = correlator.nSamples_x_chipPeriod;
end

%% tracking

[all_decSymbols, all_idTimeShift, all_idFreqShift] = demodulator.decodeOptimumShift(reader.IQsamples_float, ... 
                                                    segmentSize, shifts_delayPRN, shifts_nSamples_x_symbolPeriod, ...
                                                    shifts_nSamples_x_chipPeriod, coherenceFraction);     

figure(30)
set(gca,"ColorScale",'linear')
%N.B.: X=delay,Y=doppler
surf(demodulator.evolution(end).axis_delayPRN, ...
     demodulator.evolution(end).axis_chipPeriod, ...
     demodulator.evolution(end).trackingPeak, ...
     'EdgeColor', 'none')
% end corresponds to the latest save
%remove next lines
DelayShiftEvolution(i)=shifts_delayPRN(all_idTimeShift);
FreqShiftEvolution(i)=shifts_nSamples_x_symbolPeriod(all_idFreqShift);
i=i+1;


%% tracking splitted
inSamples = reader.IQsamples_float;

%demodulator.PRNsequence=[1 0 1 0 1 1 0 0 1 1 1 0];
%create PRN with different length due to different dopplers
PRNsampled = demodulator.getPRNFromChipPeriods(shifts_nSamples_x_chipPeriod, segmentSize);

%adapt the length of the samples to the PRN, horizontal vector
[mySamples, shifts_delayPRN] = demodulator.adaptSamplesToPRNlength(inSamples, shifts_delayPRN, size(PRNsampled, 2));

%add different delays to the PRN, sampled PRN with shifts in time and frequency
shiftedPRNsampled = demodulator.createShiftedPRN(PRNsampled, shifts_delayPRN);
%plot(shiftedPRNsampled(1:2,end-1e4:1:end)')

%in-phase & quadrature multicorrelation
corrI = demodulator.normMultiply(shiftedPRNsampled, mySamples(1, :));
corrQ = demodulator.normMultiply(shiftedPRNsampled, mySamples(2, :));

%coherent integration, over the rows there are the coherent sums
coherentCorrI = demodulator.sumOverCoherentFraction(corrI, shifts_delayPRN, ... //coherentCorrI.cols=segmentSize/coherenceFraction
                                            shifts_nSamples_x_symbolPeriod, ...
                                            coherenceFraction, segmentSize);
coherentCorrQ = demodulator.sumOverCoherentFraction(corrQ, shifts_delayPRN, ...
                                            shifts_nSamples_x_symbolPeriod, ...
                                            coherenceFraction, segmentSize);

% figure
% plot(coherentCorrI)

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
shifts_delayPRN = shifts_delayPRN - shifts_delayPRN((end+1)/2);
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
decodedSymbols(currentSymbol + 1:currentSymbol + segmentSize) = decSymbols;

%update the correlation peak with new doppler and delay estimations
new_samplesSymbolPeriod = shifts_nSamples_x_symbolPeriod(all_idFreqShift);    
advancement_startingSample = segmentSize * shifts_nSamples_x_symbolPeriod(all_idFreqShift) + ...
                             shifts_delayPRN(all_idTimeShift);
correlator.updateCorrelationPeak(new_samplesSymbolPeriod,advancement_startingSample);

currentSymbol = currentSymbol + segmentSize;
end

%%
if size(decodedSymbols, 1) > size(decodedSymbols, 2)
        decodedSymbols = decodedSymbols'; %row vector
end            
            
%minimum distance region decision, no channel inversion
decodedSymbols = int16((decodedSymbols + 1) ./ 2);

if string(num2str(decodedSymbols,'%d'))==string(num2str(packet,'%d'))
    fprintf("Message correctly received \n");
else
    fprintf("Error! Failed decoding \n");
    figure
    plot(decodedSymbols-packet,"o-")
    title("Decoding errors")
end

%%
figure
plot(1:length(DopplerCorrectionEvolution),DopplerCorrectionEvolution,'linewidth',2)
title("Doppler Correction Evolution")
axis([1 length(DopplerCorrectionEvolution) min(DopplerCorrectionEvolution) max(DopplerCorrectionEvolution)])
figure
hold on
plot(1:length(DelayShiftEvolution),DelayShiftEvolution,'linewidth',2)
title("Delay Correction Evolution")
figure
plot(1:length(FreqShiftEvolution),FreqShiftEvolution,'linewidth',2)
title("Freq Correction Evolution")

%% COMMENTI SUI RISULTATI
% Test "T_tracking_1.bin" e "T_tracking_2.bin": 
% (segmentSize=1) correlator corregge 41Hz per simbolo. 
% La DopplerCorrection decresce linearmente con il numero di simboli da
% trackare fino a quando non raggiunge la fDoppler attesa, dopo vari
% simboli (da capire perchè), si assesta attorno al valore corretto con un
% intervallo di 30Hz
