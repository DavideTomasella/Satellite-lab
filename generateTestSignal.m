%DT$ removed clearvars to hande input arguments
%clearvars
close all
addpath(".\..\")

%% signal
%mode: 0 = constant, 1 = linear, 2 = triangular, 3 = cosine, 4 = quadratic 

%TRACKING SIGNALS
PARS(1)  = struct("name","T_tracking_1",  "dstart",15.23,"dend",15.23,"deveryChip",true, "dmode",0,"envelope",0,"inNoise",0,  "envelopePhase",0,"outNoise",0,"outNoise_Length",0,"outBits_Length",0,"outPRN_Length",0,"powerReduce",0);
PARS(2)  = struct("name","T_tracking_1a", "dstart",15.23,"dend",15.23,"deveryChip",true, "dmode",0,"envelope",1,"inNoise",0,  "envelopePhase",0,"outNoise",0,"outNoise_Length",0,"outBits_Length",0,"outPRN_Length",0,"powerReduce",0);
PARS(3)  = struct("name","T_tracking_1b", "dstart",15.23,"dend",15.23,"deveryChip",true, "dmode",0,"envelope",1,"inNoise",0.1,"envelopePhase",0,"outNoise",0,"outNoise_Length",0,"outBits_Length",0,"outPRN_Length",0,"powerReduce",0);
PARS(4)  = struct("name","T_tracking_10a","dstart",15.23,"dend",18.23,"deveryChip",true, "dmode",1,"envelope",1,"inNoise",0,  "envelopePhase",0,"outNoise",0,"outNoise_Length",0,"outBits_Length",0,"outPRN_Length",0,"powerReduce",0);
PARS(5)  = struct("name","T_tracking_10b","dstart",15.23,"dend",18.23,"deveryChip",true, "dmode",1,"envelope",1,"inNoise",0.1,"envelopePhase",0,"outNoise",0,"outNoise_Length",0,"outBits_Length",0,"outPRN_Length",0,"powerReduce",0);
PARS(6)  = struct("name","T_tracking_11a","dstart",15.23,"dend",18.23,"deveryChip",true, "dmode",2,"envelope",1,"inNoise",0,  "envelopePhase",0,"outNoise",0,"outNoise_Length",0,"outBits_Length",0,"outPRN_Length",0,"powerReduce",0);
PARS(7)  = struct("name","T_tracking_11b","dstart",15.23,"dend",18.23,"deveryChip",true, "dmode",2,"envelope",1,"inNoise",0.1,"envelopePhase",0,"outNoise",0,"outNoise_Length",0,"outBits_Length",0,"outPRN_Length",0,"powerReduce",0);
PARS(8)  = struct("name","T_tracking_2",  "dstart",418.7,"dend",418.7,"deveryChip",true, "dmode",0,"envelope",0,"inNoise",0,  "envelopePhase",0,"outNoise",0,"outNoise_Length",0,"outBits_Length",0,"outPRN_Length",0,"powerReduce",0);
PARS(9)  = struct("name","T_tracking_2a", "dstart",418.7,"dend",418.7,"deveryChip",true, "dmode",0,"envelope",1,"inNoise",0,  "envelopePhase",0,"outNoise",0,"outNoise_Length",0,"outBits_Length",0,"outPRN_Length",0,"powerReduce",0);
PARS(10) = struct("name","T_tracking_2b", "dstart",418.7,"dend",418.7,"deveryChip",true, "dmode",0,"envelope",1,"inNoise",0.1,"envelopePhase",0,"outNoise",0,"outNoise_Length",0,"outBits_Length",0,"outPRN_Length",0,"powerReduce",0);
PARS(21) = struct("name","T_tracking_20", "dstart",418.7,"dend",448.7,"deveryChip",true, "dmode",1,"envelope",0,"inNoise",0,  "envelopePhase",0,"outNoise",0,"outNoise_Length",0,"outBits_Length",0,"outPRN_Length",0,"powerReduce",0);
PARS(11) = struct("name","T_tracking_20a","dstart",418.7,"dend",448.7,"deveryChip",false,"dmode",1,"envelope",1,"inNoise",0,  "envelopePhase",0,"outNoise",0,"outNoise_Length",0,"outBits_Length",0,"outPRN_Length",0,"powerReduce",0);
PARS(12) = struct("name","T_tracking_20b","dstart",418.7,"dend",448.7,"deveryChip",false,"dmode",1,"envelope",1,"inNoise",0.1,"envelopePhase",0,"outNoise",0,"outNoise_Length",0,"outBits_Length",0,"outPRN_Length",0,"powerReduce",0);
PARS(13) = struct("name","T_tracking_21a","dstart",418.7,"dend",328.7,"deveryChip",true, "dmode",1,"envelope",1,"inNoise",0,  "envelopePhase",0,"outNoise",0,"outNoise_Length",0,"outBits_Length",0,"outPRN_Length",0,"powerReduce",0);
PARS(14) = struct("name","T_tracking_21b","dstart",418.7,"dend",328.7,"deveryChip",true, "dmode",1,"envelope",1,"inNoise",0.1,"envelopePhase",0,"outNoise",0,"outNoise_Length",0,"outBits_Length",0,"outPRN_Length",0,"powerReduce",0);
PARS(15) = struct("name","T_tracking_22a","dstart",418.7,"dend",448.7,"deveryChip",true, "dmode",2,"envelope",1,"inNoise",0,  "envelopePhase",0,"outNoise",0,"outNoise_Length",0,"outBits_Length",0,"outPRN_Length",0,"powerReduce",0);
PARS(16) = struct("name","T_tracking_22b","dstart",418.7,"dend",448.7,"deveryChip",true, "dmode",2,"envelope",1,"inNoise",0.1,"envelopePhase",0,"outNoise",0,"outNoise_Length",0,"outBits_Length",0,"outPRN_Length",0,"powerReduce",0);
PARS(17) = struct("name","T_tracking_23a","dstart",418.7,"dend",448.7,"deveryChip",true, "dmode",4,"envelope",1,"inNoise",0,  "envelopePhase",0,"outNoise",0,"outNoise_Length",0,"outBits_Length",0,"outPRN_Length",0,"powerReduce",0);
PARS(18) = struct("name","T_tracking_23b","dstart",418.7,"dend",448.7,"deveryChip",true, "dmode",4,"envelope",1,"inNoise",0.1,"envelopePhase",0,"outNoise",0,"outNoise_Length",0,"outBits_Length",0,"outPRN_Length",0,"powerReduce",0);
PARS(19) = struct("name","T_tracking_24a","dstart",418.7,"dend",448.7,"deveryChip",true, "dmode",3,"envelope",1,"inNoise",0,  "envelopePhase",0,"outNoise",0,"outNoise_Length",0,"outBits_Length",0,"outPRN_Length",0,"powerReduce",0);
PARS(20) = struct("name","T_tracking_24b","dstart",418.7,"dend",448.7,"deveryChip",true, "dmode",3,"envelope",1,"inNoise",0.1,"envelopePhase",0,"outNoise",0,"outNoise_Length",0,"outBits_Length",0,"outPRN_Length",0,"powerReduce",0);

%AQUISITION SIGNALS
PARS(21)  = struct("name","T_acquisition_1", "dstart",15.23,"dend",18.23,"deveryChip",true, "dmode",1,"envelope",1,"inNoise",0,  "envelopePhase",0,  "outNoise",0,  "outNoise_Length",0,  "outBits_Length",0, "outPRN_Length",0,  "powerReduce",0);
PARS(22)  = struct("name","T_acquisition_1a","dstart",15.23,"dend",18.23,"deveryChip",true, "dmode",1,"envelope",1,"inNoise",0,  "envelopePhase",0.5,"outNoise",0,  "outNoise_Length",0,  "outBits_Length",0, "outPRN_Length",0,  "powerReduce",0);
PARS(23)  = struct("name","T_acquisition_1b","dstart",15.23,"dend",18.23,"deveryChip",true, "dmode",1,"envelope",1,"inNoise",0.8,"envelopePhase",0,  "outNoise",0,  "outNoise_Length",0,  "outBits_Length",0, "outPRN_Length",0,  "powerReduce",0);
PARS(24)  = struct("name","T_acquisition_2", "dstart",15.23,"dend",18.23,"deveryChip",true, "dmode",1,"envelope",1,"inNoise",0,  "envelopePhase",0,  "outNoise",0.4,"outNoise_Length",100,"outBits_Length",0, "outPRN_Length",0,  "powerReduce",0);
PARS(25)  = struct("name","T_acquisition_2a","dstart",15.23,"dend",18.23,"deveryChip",true, "dmode",1,"envelope",1,"inNoise",0.4,"envelopePhase",0,  "outNoise",0.4,"outNoise_Length",100,"outBits_Length",0, "outPRN_Length",0,  "powerReduce",0);
PARS(26)  = struct("name","T_acquisition_2b","dstart",15.23,"dend",18.23,"deveryChip",true, "dmode",1,"envelope",1,"inNoise",0.8,"envelopePhase",0,  "outNoise",0.4,"outNoise_Length",100,"outBits_Length",0, "outPRN_Length",0,  "powerReduce",0);
PARS(27)  = struct("name","T_acquisition_2c","dstart",15.23,"dend",18.23,"deveryChip",true, "dmode",1,"envelope",1,"inNoise",0.8,"envelopePhase",0,  "outNoise",0.4,"outNoise_Length",100,"outBits_Length",0, "outPRN_Length",0,  "powerReduce",1);
PARS(28)  = struct("name","T_acquisition_2d","dstart",15.23,"dend",18.23,"deveryChip",true, "dmode",1,"envelope",1,"inNoise",0.8,"envelopePhase",0,  "outNoise",0.4,"outNoise_Length",100,"outBits_Length",0, "outPRN_Length",0,  "powerReduce",2);
PARS(29)  = struct("name","T_acquisition_2e","dstart",15.23,"dend",18.23,"deveryChip",true, "dmode",1,"envelope",1,"inNoise",0.8,"envelopePhase",0,  "outNoise",0.4,"outNoise_Length",100,"outBits_Length",0, "outPRN_Length",0,  "powerReduce",3);
PARS(30)  = struct("name","T_acquisition_2f","dstart",15.23,"dend",18.23,"deveryChip",true, "dmode",1,"envelope",1,"inNoise",0.8,"envelopePhase",0,  "outNoise",0.4,"outNoise_Length",100,"outBits_Length",0, "outPRN_Length",0,  "powerReduce",4);
PARS(31)  = struct("name","T_acquisition_2g","dstart",15.23,"dend",18.23,"deveryChip",true, "dmode",1,"envelope",1,"inNoise",0.8,"envelopePhase",0,  "outNoise",0.4,"outNoise_Length",100,"outBits_Length",0, "outPRN_Length",0,  "powerReduce",5);
PARS(32)  = struct("name","T_acquisition_3", "dstart",15.23,"dend",18.23,"deveryChip",true, "dmode",1,"envelope",1,"inNoise",0.8,"envelopePhase",0,  "outNoise",0.8,"outNoise_Length",250,"outBits_Length",0, "outPRN_Length",0,  "powerReduce",0);
PARS(33)  = struct("name","T_acquisition_3a","dstart",15.23,"dend",18.23,"deveryChip",true, "dmode",1,"envelope",1,"inNoise",0,  "envelopePhase",0,  "outNoise",0.8,"outNoise_Length",250,"outBits_Length",0, "outPRN_Length",100,"powerReduce",0);
PARS(34)  = struct("name","T_acquisition_3b","dstart",15.23,"dend",18.23,"deveryChip",true, "dmode",1,"envelope",1,"inNoise",0,  "envelopePhase",0,  "outNoise",0.8,"outNoise_Length",250,"outBits_Length",3, "outPRN_Length",400,"powerReduce",0);
PARS(35)  = struct("name","T_acquisition_3c","dstart",15.23,"dend",18.23,"deveryChip",true, "dmode",1,"envelope",1,"inNoise",0.8,"envelopePhase",0,  "outNoise",0.8,"outNoise_Length",250,"outBits_Length",13,"outPRN_Length",400,"powerReduce",0);
PARS(36)  = struct("name","T_acquisition_3d","dstart",15.23,"dend",18.23,"deveryChip",true, "dmode",1,"envelope",1,"inNoise",0.8,"envelopePhase",0,  "outNoise",0.8,"outNoise_Length",250,"outBits_Length",18,"outPRN_Length",400,"powerReduce",0);

PARS(37)  = struct("name","T_acquisition_test","dstart",15.23,"dend",18.23,"deveryChip",true, "dmode",1,"envelope",1,"inNoise",0,  "envelopePhase",0.5,"outNoise",0,  "outNoise_Length",0,  "outBits_Length",0, "outPRN_Length",0,  "powerReduce",0);
PARS(38)  = struct("name","T_tracking_1c", "dstart",15.23,"dend",15.23,"deveryChip",true, "dmode",0,"envelope",1,"inNoise",0.1,"envelopePhase",0,"outNoise",0,"outNoise_Length",0,   "outBits_Length",0,"outPRN_Length",0, "powerReduce",1);
PARS(39)  = struct("name","T_tracking_1d", "dstart",15.23,"dend",15.53,"deveryChip",true, "dmode",1,"envelope",1,"inNoise",1,  "envelopePhase",0,"outNoise",1,"outNoise_Length",1203,"outBits_Length",0,"outPRN_Length",45,"powerReduce",6);
PARS(40)  = struct("name","T_tracking_1e", "dstart",15.23,"dend",15.53,"deveryChip",true, "dmode",1,"envelope",1,"inNoise",0,"envelopePhase",0,"outNoise",1,"outNoise_Length",1203,"outBits_Length",0,"outPRN_Length",45,"powerReduce",13);

p = 19;

if ~exist('DEBUG',"var")
    DEBUG = false;
end
if exist('TESTPARS',"var")
    PARS(p) = TESTPARS;
end
if ~exist('signalDirectory',"var")
    signalDirectory = "binData/testSignals";
end

%for p=1:36

%% Configuration settings and files
inout = InOutInterface();
% directories
indir = "inData";
inout.configCreateSettings(indir);
outdir = "outData";
inout.configSaveResults(outdir);

if p == 23
    filename = "in_PRN2.json";
else
    filename = "in0.json";
end
PRNfilename = "PRNpattern.json";
settings = inout.createSettings(filename,PRNfilename);

reader = SignalManager();
reader.configReadFile(signalDirectory, "nine.bin", inout.settings.quantizationBits);
date = datestr(now, '_yymmdd_HHMMSS');
%outputFileName = strcat("T_tracking_1", date, ".bin");
outputFileName = strcat(PARS(p).name,".bin");
 

%% Create bits (message)
%Corrected pkt
if p == 37
    packet = repmat([0  1], 1, 200);
else
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
end
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
BITS = repelem(packet, 1, inout.settings.nPRN_x_Symbol);

%% Add random bits in front and tail
preBLength = PARS(p).outBits_Length;
postBLength = PARS(p).outBits_Length;
BITS = [randi([0,1], 1, preBLength) BITS randi([0,1], 1, postBLength)];

%% Create symbols
SYMBOLS = (2 * BITS - 1);
%DT$ add always 1 symbol at the end (TEMP FIX FOR INTERP1 SKIPPING LAST SYMBOL)
SYMBOLS = [SYMBOLS 1];

%% Multiply per PRN
PRNsequence = 2 * inout.PRNcode - 1; 
SEQUENCE = reshape(SYMBOLS.*PRNsequence',1,[]);

%% Add PRN peaces
prePLength = PARS(p).outPRN_Length; %max 4092
postPLength = PARS(p).outPRN_Length;
SEQUENCE = [PRNsequence(1:prePLength) SEQUENCE PRNsequence(1:postPLength)];

%% Creation of base signal (strecthed due to doppler)
 % array 1 campione per chip -> doppler variabile = lineare
               % aggiungere rumore al vettore
               % trasformarlo anche questo con upsampling  ("linear")per 
               % usarlo in exp(1i*2*pi*fdoppler*..)
%SEQUENCE=[1 0 1 0 1 1 0 0 1 1 1 0];

advancedeveryChip = PARS(p).deveryChip;
mode = PARS(p).dmode; 
%mode=0 : constant doppler
%mode=1 : linear increase
%mode=2 : triangolar increase
%mode=3 : cosine increase
%mode=4 : quadratic increase
dopStart = PARS(p).dstart;
dopEnd = PARS(p).dend;
%create valueDoppler from specification
nDopplerSteps = length(SEQUENCE);
if ~advancedeveryChip %CHANGE DOPPLER EVERY CHIP
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

if ~advancedeveryChip
    valueDoppler = repelem(valueDoppler,1,length(inout.PRNcode)*inout.settings.nPRN_x_Symbol);
end
newChipRate = inout.settings.chipRate + valueDoppler; %add doppler

%Time+doppler vector creation
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
outSigma = PARS(p).outNoise;    % noise variance 
preSLength =PARS(p).outNoise_Length; %samples
postSLength = PARS(p).outNoise_Length;
genSIGNAL = [outSigma * randn(preSLength, 1); genSIGNAL; outSigma * randn(postSLength, 1)];
genDOPPLER = [genDOPPLER(1) * ones(preSLength, 1); genDOPPLER; genDOPPLER(end) * ones(preSLength, 1)];
t = (0:1:(length(t) + preSLength + postSLength - 1))' / inout.settings.fSampling;

%% Set signal power
maxAmplitude = 2 ^ (inout.settings.quantizationBits - 2); % gain of the signal, otherwise the in16 matrix results made of 1,0 and -1
attenuation = PARS(p).powerReduce;
genSIGNAL = maxAmplitude / (2 ^ attenuation) * genSIGNAL;

%% Add noise
inSigma = PARS(p).inNoise;    % noise variance 
noise = inSigma / sqrt(2) * randn(1, length(t))';
genSIGNAL = genSIGNAL + maxAmplitude * noise;

%% Add doopler envelope and orthogonal noise
phase = PARS(p).envelopePhase;
env = PARS(p).envelope; %If 1 add envelope, 0 just noise
genSIGNAL = genSIGNAL.*exp(env*1i*(2*pi*genDOPPLER.*t+phase));
orthoNOISE = 1i*inSigma / sqrt(2) * randn(1, length(t))'.*exp(env*1i*(2*pi*genDOPPLER.*t+phase));
genSIGNAL = genSIGNAL + maxAmplitude * orthoNOISE;
reader.IQsamples_float = [real(genSIGNAL) imag(genSIGNAL)];

%% Plot
if DEBUG
    figure(60)
    plot(reader.IQsamples(:,1))
    xlim([280e3 281e3])
    %figure(61)
    %plot(reader.IQsamples(:,1))
    %xlim([1 1e5])
    %figure(62)
    %plot(reader.IQsamples(:,1))
    figure(63)
    plot(genDOPPLER)
    pause(1)
end
%PLOTTING
% h = figure(1);
% plot(t,genDOPPLER);
% xlabel("Time[s]");
% ylabel("Frequency[Hz]");
% title("Sinsusoidal Doppler Variation");
% savePdf(h,"sinusoid");
%% Save binary file
%reader.saveToBynaryFile(reader.IQsamples,outputFileName);
sprintf("Generation test signal completed.")

%end