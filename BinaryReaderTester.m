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
reader.configReadFile("inData","nine.bin",inout.settings.quantizationBits);
%EXAMPLE this line reads a window of 10 samples starting from
%the 5th sample(included). The samples are stored in the field IQsamples of
%the reader object
reader.readFile(5,10);
samples = reader.IQsamples;
%EXAMPLE this lines these lines create a file called prova.bin in the
%current folder, if already present appends the contend specified in the
%column vector samples to the file
status = reader.appendToBynaryFile("prova.bin",samples);
status1 = reader.appendToBynaryFile("prova.bin",samples);
%New binary reader for testing
reader1 = BinaryReader();
reader1.configReadFile(".","prova.bin",inout.settings.quantizationBits);
reader1.readFile(1,20);
samples1 = reader1.IQsamples;