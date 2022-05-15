clearvars
close all

%% Configuration settings and files
inout = InOutInterface();
% directories
indir = "inData";
inout.configCreateSettings(indir);
outdir = "outData";
inout.configSaveResults(outdir);

filename = "in0.json";
PRNfilename = "PRNpattern.json";
settings = inout.createSettings(filename,PRNfilename);

reader = BinaryReader();
reader.configReadFile("binData/testSignals", "nine.bin", inout.settings.quantizationBits);
date = datestr(now, '_yymmdd_HHMMSS');
outputFileName = strcat("T_tracking_4", date, ".bin");
%outputFileName = strcat("signal_test_A.bin");
 

%% Create symbols
% packet = [0  1  0  1  1  0  0  0  ...
%           0  0  0  0  0  0  0  1  ...
%           0  0  1  1  0  0  0  0  ...
%           0  0  0  0  0  0  0  0  ...
%           0  0  0  0  1  1  1  1  ...
%           1  1  1  1  1  1  1  1  ...
%           1  1  1  0  0  0  0  1  ...
%           1  0  1  1  0  1  0  1  ...
%           1  0  0  0  0  1  0  1  ...
%           1  1  0  0  0  0  0  0];
packet = [0  1  0  1  1  0  0  0  ...
          0  0  0  0  0  0  0  1  ...
          0  0  1  1];
BITS = repelem(packet,1,inout.settings.nPRN_x_Symbol);
SYMBOLS = (2 * BITS - 1);

%% Add random symbols in front and tail
preMLength = 0;
postMLength = 0;
SYMBOLS = [rand(1,preMLength), SYMBOLS, rand(1,postMLength)];

%% Multiply per PRN
SEQUENCE = reshape(SYMBOLS.*inout.PRNcode',1,[]);

%% Creation of base signal (strecthed due to doppler)
 % array 1 campione per chip -> doppler variabile = lineare
               % aggiungere rumore al vettore
               % trasformarlo anche questo con upsampling  ("linear")per 
               % usarlo in exp(1i*2*pi*fdoppler*..)
addLinearChirp = true;
if addLinearChirp 
    dopMin = 418.7;
    dopMax = 448.7;
    stepDoppler = (dopMax-dopMin)/(length(SYMBOLS)-1);
    advanceEveryChip = true;
    if advanceEveryChip
        stepDoppler = (dopMax-dopMin)/(length(SEQUENCE)-1);
        fd = (0:length(SEQUENCE)-1)*stepDoppler;
        newChipRate = inout.settings.chipRate + fd; %add doppler
    else
        stepDoppler = (dopMax-dopMin)/(length(SYMBOLS)-1);
        fd = (0:length(SYMBOLS)-1)*stepDoppler;
        fdd = repelem(fd,1,length(inout.PRNcode)*inout.settings.nPRN_x_Symbol);
        newChipRate = inout.settings.chipRate + fdd; %add doppler
    end
else
    dop_const = 418.7;
    dop_const_vect = ones(1,length(SEQUENCE))*dop_const;
    newChipRate = inout.settings.chipRate + dop_const_vect; %add doppler
end
tchip = ((0:length(SEQUENCE)-1)./ newChipRate)';
fdoppler = (newChipRate-inout.settings.chipRate)';
%IMPORTANT: verificare che siano uguali i due vettori di tempi
%upsampling = inout.settings.fSampling / newChipRate;
%t = (0:upsampling*length(SEQUENCE)-1)' / inout.settings.fSampling;
t = (0:1/inout.settings.fSampling:tchip(end))';
genSIGNAL = interp1(tchip,SEQUENCE,t,"previous");
genDoppler = interp1(tchip,fdoppler,t,"previous");
%% Add front and tail samples
sigmaPP = 0.0;    % noise variance 
preSLength = 0; %samples
postSLength = 0;
genSIGNAL = [sigmaPP * randn(1,preMLength), genSIGNAL, sigmaPP * randn(1,postMLength)];

%% Add noise
sigma = 0;    % noise variance 
noise = sigma / sqrt(2) * randn(1, length(t))';
genSIGNAL = genSIGNAL + noise;

%% Set signal power
amplitude = 2 ^ (inout.settings.quantizationBits - 2); % gain of the signal, otherwise the in16 matrix results made of 1,0 and -1
genSIGNAL = amplitude * genSIGNAL;

%% Add doopler envelope and orthogonal noise
delay = 0;
env = 0; %If 1 add envelope, 0 just noise
genSIGNAL = genSIGNAL.*exp(env*1i*2*pi*genDoppler.*(t-delay));
orthoNOISE = 1i*sigma / sqrt(2) * randn(1, length(t))'.*exp(env*1i*2*pi*genDoppler.*(t-delay));
genSIGNAL = genSIGNAL + orthoNOISE;
reader.IQsamples_float = [real(genSIGNAL) imag(genSIGNAL)];

%% Plot
figure(10)
plot(reader.IQsamples(:,1))
xlim([1 1e3])
figure(11)
plot(reader.IQsamples(:,1))
xlim([1 1e5])

%% Save binary file
reader.saveToBynaryFile(reader.IQsamples,outputFileName);
