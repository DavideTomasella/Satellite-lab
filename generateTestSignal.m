clearvars
close all

%% signal
%mode: 0 = constant, 1 = linear, 2 = triangular, 3 = cosine, 4 = quadratic 
clear vars
PARS(1) = struct("name","T_tracking_1","dstart",15.23,"dend",15.23,"everyChip",true,"mode",0,"envelope",0,"noise",0);
PARS(2) = struct("name","T_tracking_1a","dstart",15.23,"dend",15.23,"everyChip",true,"mode",0,"envelope",1,"noise",0);
PARS(3) = struct("name","T_tracking_1b","dstart",15.23,"dend",15.23,"everyChip",true,"mode",0,"envelope",1,"noise",1);
PARS(4) = struct("name","T_tracking_10a","dstart",15.23,"dend",18.23,"everyChip",true,"mode",1,"envelope",1,"noise",0);
PARS(5) = struct("name","T_tracking_10b","dstart",15.23,"dend",18.23,"everyChip",true,"mode",1,"envelope",1,"noise",1);
PARS(6) = struct("name","T_tracking_11a","dstart",15.23,"dend",18.23,"everyChip",true,"mode",2,"envelope",1,"noise",0);
PARS(7) = struct("name","T_tracking_11b","dstart",15.23,"dend",18.23,"everyChip",true,"mode",2,"envelope",1,"noise",1);
PARS(8) = struct("name","T_tracking_2","dstart",418.7,"dend",418.7,"everyChip",true,"mode",0,"envelope",0,"noise",0);
PARS(9) = struct("name","T_tracking_2a","dstart",418.7,"dend",418.7,"everyChip",true,"mode",0,"envelope",1,"noise",0);
PARS(10) = struct("name","T_tracking_2b","dstart",418.7,"dend",418.7,"everyChip",true,"mode",0,"envelope",1,"noise",1);
PARS(21) = struct("name","T_tracking_20","dstart",418.7,"dend",448.7,"everyChip",true,"mode",1,"envelope",0,"noise",0);
PARS(11) = struct("name","T_tracking_20a","dstart",418.7,"dend",448.7,"everyChip",false,"mode",1,"envelope",1,"noise",0);
PARS(12) = struct("name","T_tracking_20b","dstart",418.7,"dend",448.7,"everyChip",false,"mode",1,"envelope",1,"noise",1);
PARS(13) = struct("name","T_tracking_21a","dstart",418.7,"dend",328.7,"everyChip",true,"mode",1,"envelope",1,"noise",0);
PARS(14) = struct("name","T_tracking_21b","dstart",418.7,"dend",328.7,"everyChip",true,"mode",1,"envelope",1,"noise",1);
PARS(15) = struct("name","T_tracking_22a","dstart",418.7,"dend",448.7,"everyChip",true,"mode",2,"envelope",1,"noise",0);
PARS(16) = struct("name","T_tracking_22b","dstart",418.7,"dend",448.7,"everyChip",true,"mode",2,"envelope",1,"noise",1);
PARS(17) = struct("name","T_tracking_23a","dstart",418.7,"dend",448.7,"everyChip",true,"mode",4,"envelope",1,"noise",0);
PARS(18) = struct("name","T_tracking_23b","dstart",418.7,"dend",448.7,"everyChip",true,"mode",4,"envelope",1,"noise",1);
PARS(19) = struct("name","T_tracking_24a","dstart",418.7,"dend",448.7,"everyChip",true,"mode",3,"envelope",1,"noise",0);
PARS(20) = struct("name","T_tracking_24b","dstart",418.7,"dend",448.7,"everyChip",true,"mode",3,"envelope",1,"noise",1);




p = 13;
PARS(p);
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
%outputFileName = strcat("T_tracking_1", date, ".bin");
outputFileName = strcat(PARS(p).name,".bin");
 

%% Create symbols
%Corrected pkt
packet = [0  1  0  1  1  0  0  0  ...
            0  0  0  0  0  0  0  1  ...
            0  0  1  1  0  1  0  1  ...
            0  0  1  1  0  0  0  1  ...
            1  1  0  0  1  1  0  0  ...
            1  0  1  1  1  0  1  1  ...
            1  0  0  0  0  1  1  0  ...
            1  1  1  1  0  0  1  0  ...
            0  0  1  0  1  0  0  0  ...
            0  1  0  0  0  0  0  0];
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
%packet = [0  1  0  1  1  0  0  0  ...
%          0  0  0  0  0  0  0  1  ...
%          0  0  1  1];
BITS = repelem(packet,1,inout.settings.nPRN_x_Symbol);
SYMBOLS = (2 * BITS - 1);

%% Add random symbols in front and tail
preMLength = 0;
postMLength = 0;
%DT$ add always 1 symbol at the end (TEMP FIX FOR INTERP1 SKIPPING LAST SYMBOL)
SYMBOLS = [rand(1,preMLength), SYMBOLS, rand(1,postMLength) 1];

%% Multiply per PRN
PRNsequence = 2 * inout.PRNcode - 1; 
SEQUENCE = reshape(SYMBOLS.*PRNsequence',1,[]);

%% Creation of base signal (strecthed due to doppler)
 % array 1 campione per chip -> doppler variabile = lineare
               % aggiungere rumore al vettore
               % trasformarlo anche questo con upsampling  ("linear")per 
               % usarlo in exp(1i*2*pi*fdoppler*..)
%SEQUENCE=[1 0 1 0 1 1 0 0 1 1 1 0];

advanceEveryChip = PARS(p).everyChip;
mode = PARS(p).mode; 
%mode=0 : constant doppler
%mode=1 : linear increase
%mode=2 : triangolar increase
%mode=3 : cosine increase
%mode=4 : quadratic increase
dopStart = PARS(p).dstart;
dopEnd = PARS(p).dend;
%create valueDoppler from specification
nDopplerSteps = length(SEQUENCE);
if ~advanceEveryChip %CHANGE DOPPLER EVERY CHIP
    nDopplerSteps = nDopplerSteps / (length(inout.PRNcode) * inout.settings.nPRN_x_Symbol);
end
timeFraction = (0:nDopplerSteps-1) / nDopplerSteps;
if mode == 0 %constant
    value = dopStart;
    valueDoppler = 0 * timeFraction + value;
elseif mode == 1 %linear
    value = (dopEnd-dopStart);
    valueDoppler = timeFraction * value + dopStart;
elseif mode == 2 %triangular
    value = (dopEnd-dopStart);
    valueDoppler = 2 * timeFraction .* (timeFraction < 0.5) * value + ...
                   2 * (1 - timeFraction) .* (timeFraction >= 0.5) * value + ... 
                   dopStart;
    % since time is 0 -> 1, to create the square we just .^2 and normalize
elseif mode == 3 %cosine
    timeFraction = -pi + timeFraction * 2 * pi;
    value = (dopEnd-dopStart);
    valueDoppler = 0.5 * (1 + cos(timeFraction)) * value + ...
                   dopStart;
elseif mode == 4 %quadratic
    value = (dopEnd-dopStart);
    valueDoppler = timeFraction .^2 * value + dopStart;
else
    value = dopStart;
    valueDoppler = 0 * timeFraction + value;
end

%Triangular doppler
%           stepDoppler = 2*((dopMax-dopMin)/(length(SEQUENCE)-1));
%           fd1 = (0:round((length(SEQUENCE)-1)/2))*stepDoppler;
%           fd2 = flip(fd1(1:(length(SEQUENCE)-length(fd1))));
%           fd = [fd1 fd2];

%construct the doppler frequencies and chipRate vectors from valueDoppler

if ~advanceEveryChip
    valueDoppler = repelem(valueDoppler,1,length(inout.PRNcode)*inout.settings.nPRN_x_Symbol);
end
newChipRate = inout.settings.chipRate + valueDoppler; %add doppler

tchip = ((0:length(SEQUENCE)-1)./ newChipRate)';
fdoppler = (newChipRate-inout.settings.chipRate)';
%IMPORTANT: verificare che siano uguali i due vettori di tempi
%upsampling = inout.settings.fSampling / newChipRate;
%t = (0:upsampling*length(SEQUENCE)-1)' / inout.settings.fSampling;
t = (0:1/inout.settings.fSampling:tchip(end))';
genSIGNAL = interp1(tchip,SEQUENCE,t,"previous");
genDOPPLER = interp1(tchip,fdoppler,t,"previous");
%pspectrum(genSIGNAL,"persistence","FrequencyLimits",[0 0.5],"")
%eyediagram(genSIGNAL(1:40920),int16(1/inout.settings.chipRate*inout.settings.fSampling))

%% Add front and tail samples
sigmaPP = 0.0;    % noise variance 
preSLength = 0; %samples
postSLength = 0;
genSIGNAL = [sigmaPP * randn(1,preMLength), genSIGNAL, sigmaPP * randn(1,postMLength)];

%% Add noise
sigma = PARS(p).noise;    % noise variance 
noise = sigma / sqrt(2) * randn(1, length(t))';
genSIGNAL = genSIGNAL + noise;

%% Set signal power
amplitude = 2 ^ (inout.settings.quantizationBits - 2); % gain of the signal, otherwise the in16 matrix results made of 1,0 and -1
genSIGNAL = amplitude * genSIGNAL;

%% Add doopler envelope and orthogonal noise
delay = 0;
phase = 0;   % test for a phase
env = PARS(p).envelope; %If 1 add envelope, 0 just noise
genSIGNAL = genSIGNAL.*exp(env*1i*2*pi*genDOPPLER.*(t-delay)+phase);
orthoNOISE = 1i*sigma / sqrt(2) * randn(1, length(t))'.*exp(env*1i*2*pi*genDOPPLER.*(t-delay));
genSIGNAL = genSIGNAL + orthoNOISE;
reader.IQsamples_float = [real(genSIGNAL) imag(genSIGNAL)];

%% Plot
figure(10)
plot(reader.IQsamples(:,1))
xlim([1 1e3])
figure(11)
plot(reader.IQsamples(:,1))
xlim([1 1e5])
figure(12)
plot(reader.IQsamples(:,1))
figure(13)
plot(genDOPPLER)

%% Save binary file
reader.saveToBynaryFile(reader.IQsamples,outputFileName);
