%%analyzeMessage tester%%
%Settings/config
clearvars; close all; clc;

inout = InOutInterface();
%Define in-out directories
inout.configCreateSettings("inData");
inout.configSaveResults("outData");
%Read settings
inout.createSettings("in0.json");


dem = Demodulator();

dem.configMessageAnalyzer(inout.settings.CRCpolynomial,inout.settings.SV_PRN_ID, ...
                          inout.settings.SYNCpattern, inout.settings.SVIDlength, ...
                          inout.settings.MIDlength,inout.settings.MBODYlength_ACK, ...
                          inout.settings.CRClength, inout.settings.TAILpattern);

%% 
%Parameters%

SV_ID=[0 0 0 1 0 1]; %correct SV_ID should be [0 0 0 0 0 1] N.B.: 6 elements
%SyncPattern=inout.settings.SYNCpattern-'0'; %correct sync N.B.: 10
%elements
SyncPattern=[0 1 0 1 1 0 0 0 0 1]; 
M_ID=[0 0 1 1 0]; 
M_body=[1, zeros(1,15), ones(1,14)];
TailPattern=[0 0 0 0 0 1];

CRCpoly=[1 0 0 0 0 0 1 0 0 0 1 1 1 0 1 1 1 0 1 0 1 0 0 1]; %not ACKed
%CRCpoly=dem.calcChecksum([M_ID M_body zeros(1,inout.settings.CRClength)]);

%% 
%Message generation (symbols)%
message=[SyncPattern, SV_ID, M_ID, M_body, CRCpoly, TailPattern];
correctMessage=[0  1  0  1  1  0  0  0  ...
                0  0  0  0  0  0  0  1  ...
                0  0  1  1  0  0  0  0  ...
                0  0  0  0  0  0  0  0  ...
                0  0  0  0  1  1  1  1  ...
                1  1  1  1  1  1  1  1  ...
                1  1  1  0  0  0  0  1  ...
                1  0  1  1  0  1  0  1  ...
                1  0  0  0  0  1  0  1  ...
                1  1  0  0  0  0  0  0];
decodedSymbols=message.*2-1;
%decodedSymbols=correctMessage.*2-1;

outData=inout;
dem.analyzeMessage(decodedSymbols, outData);
outData.results
%saved = inout.saveResults("test_results.json");