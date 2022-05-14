%CODE TO CREATE A BINARY FILE called nine.bin in the current folder
% fileID = fopen('nine.bin','w');
% a = 1:1:30;
% a_16= int16(a);
% fwrite(fileID,a,'int16');
% fclose(fileID);

%Copy paste spudorato

%RECEIVER MAIN FILE (CONFIGURATION AND PROCESS ROUTINE)
clearvars
close all
addpath(".\..\")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: SETTINGS ACQUITION        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inout = InOutInterface();
%Define in-out directories
inout.configCreateSettings("inData");
inout.configSaveResults("outData");
%Read settings
inout.createSettings("in0.json");

%creating the reader object
reader = BinaryReader();
%The test binary file is ./inData/nine.bin
reader.configReadFile("binData","noise_doppler_sinc_pattern.bin",inout.settings.quantizationBits);
%EXAMPLE this line reads a window of 6 samples starting from
%the 1th sample(included) (skip=0 samples). The samples are stored in the field IQsamples of
%the reader object
nS = reader.nSamples;

% samples = reader.IQsamples;
% %EXAMPLE this lines these lines create a file called prova.bin in the
% %current folder, if already present appends the contend specified in the
% %column vector samples to the file
% reader.saveToBynaryFile(samples,"prova.bin",false);
% reader.saveToBynaryFile(samples,"prova.bin",true);
% %New binary reader for testing
% reader1 = BinaryReader();
% reader1.configReadFile("binData","signal_test1.bin",inout.settings.quantizationBits);
% reader1.readFile(0,12);
% samples1 = reader1.IQsamples;
% if sum(samples1-[samples;samples],"all")>0 disp("Error in saving or reading binary file"); end