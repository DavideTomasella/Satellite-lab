close all; clc;
PRNsequence=repelem([1 -1],10);
e_nSamples_x_chipPeriod=10.5;

PRNinterp = interp1(1:length(PRNsequence), PRNsequence, ...
    1:1 / e_nSamples_x_chipPeriod:length(PRNsequence), "previous"); %upsampling 

figure
hold on
stem(1:1 / e_nSamples_x_chipPeriod:length(PRNsequence),PRNinterp,'r')
stem(1:length(PRNsequence),PRNsequence,'b')
hold off