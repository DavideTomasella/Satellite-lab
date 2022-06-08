addpath("reports\")
%REPORT=load("reports\resultsNONOISE\REPORT_backup3.mat");
REPORT=load("reports\resultsNOISE\REPORT_K_backup1010.mat");
clear sumMatrix
clear maxMatrix
lAtten = 10;
lDopp = 10;
sumMatrix = zeros(lAtten,lDopp,"int32");
bestMatrix = zeros(lAtten,lDopp,"int32");
idBestMatrix = string(bestMatrix);
nBestMatrix = zeros(lAtten,lDopp,"int32");
fields = string(fieldnames(REPORT)');

%output matrices and relative filters
nofilterMatrix = zeros(lAtten,lDopp,"int32");
nofilter = @(x) contains(x,'f0');
filter1Matrix = zeros(lAtten,lDopp,"int32");
filter1 = @(x) contains(x,'f1');
filter3Matrix = zeros(lAtten,lDopp,"int32");
filter3 = @(x) contains(x,'f3');
resol1Matrix = zeros(lAtten,lDopp,"int32");
resol1 = @(x) contains(x,'d100s');
resol5Matrix = zeros(lAtten,lDopp,"int32");
resol5 = @(x) contains(x,'d500s');
resol10Matrix = zeros(lAtten,lDopp,"int32");
resol10 = @(x) contains(x,'d1000');
resol20Matrix = zeros(lAtten,lDopp,"int32");
resol20 = @(x) contains(x,'d2000');
resol50Matrix = zeros(lAtten,lDopp,"int32");
resol50 = @(x) contains(x,'d5000');
segm1Matrix = zeros(lAtten,lDopp,"int32");
segm1 = @(x) contains(x,'s1c');
segm2Matrix = zeros(lAtten,lDopp,"int32");
segm2 = @(x) contains(x,'s2');
segm5Matrix = zeros(lAtten,lDopp,"int32");
segm5 = @(x) contains(x,'s5');
segm10Matrix = zeros(lAtten,lDopp,"int32");
segm10 = @(x) contains(x,'s10');
coher1Matrix = zeros(lAtten,lDopp,"int32");
coher1 = @(x) contains(x,'c1');
coher3Matrix = zeros(lAtten,lDopp,"int32");
coher3 = @(x) contains(x,'c3');
coher5Matrix = zeros(lAtten,lDopp,"int32");
coher5 = @(x) contains(x,'c5');

%exclude part of the data
excludeData = @(f,fun) f(~fun(f));
exclude = true;
namePdf = "tmp";
if exclude
    fields = excludeData(fields,nofilter);
    fields = excludeData(fields,filter3);
    %fields = excludeData(fields,filter1);
    fields = excludeData(fields,resol50);
    fields = excludeData(fields,resol20);
    fields = excludeData(fields,resol10);
    fields = excludeData(fields,resol5);
    %fields = excludeData(fields,resol1);
    fields = excludeData(fields,segm10);
    fields = excludeData(fields,segm5);
    %fields = excludeData(fields,segm2);
    fields = excludeData(fields,segm1);
    fields = excludeData(fields,coher5);
    fields = excludeData(fields,coher3);
    %fields = excludeData(fields,coher1);
end

for sec = fields
    sumMatrix = sumSub(sumMatrix, int32(REPORT.(sec)));
    [bestMatrix, idBestMatrix, nBestMatrix] = ...
        maxSub(bestMatrix, idBestMatrix, nBestMatrix, int32(REPORT.(sec)), sec);
    if nofilter(sec) nofilterMatrix = sumSub(nofilterMatrix,int32(REPORT.(sec))); end
    if filter1(sec) filter1Matrix = sumSub(filter1Matrix, int32(REPORT.(sec))); end
    if filter3(sec) filter3Matrix = sumSub(filter3Matrix, int32(REPORT.(sec))); end
    if resol1(sec) resol1Matrix = sumSub(resol1Matrix, int32(REPORT.(sec))); end
    if resol5(sec) resol5Matrix = sumSub(resol5Matrix, int32(REPORT.(sec))); end
    if resol10(sec) resol10Matrix = sumSub(resol10Matrix, int32(REPORT.(sec))); end
    if resol20(sec) resol20Matrix = sumSub(resol20Matrix, int32(REPORT.(sec))); end
    if resol50(sec) resol50Matrix = sumSub(resol50Matrix, int32(REPORT.(sec))); end
    if segm1(sec) segm1Matrix = sumSub(segm1Matrix, int32(REPORT.(sec))); end
    if segm2(sec) segm2Matrix = sumSub(segm2Matrix, int32(REPORT.(sec))); end
    if segm5(sec) segm5Matrix = sumSub(segm5Matrix, int32(REPORT.(sec))); end
    if segm10(sec) segm10Matrix = sumSub(segm10Matrix, int32(REPORT.(sec))); end
    if coher1(sec) coher1Matrix = sumSub(coher1Matrix, int32(REPORT.(sec))); end
    if coher3(sec) coher3Matrix = sumSub(coher3Matrix, int32(REPORT.(sec))); end
    if coher5(sec) coher5Matrix = sumSub(coher5Matrix, int32(REPORT.(sec))); end
end
disp(flipud(idBestMatrix))
h=figure(8);
movegui("southwest")
setFigure(sumMatrix/length(fields),"ALL DATA")
%savePdf(h,"reports/" + namePdf,true)
figure(9)
movegui("south")
setFigure(bestMatrix,"BEST DATA")
figure(10)
movegui("southeast")
setFigure(nBestMatrix,"SUCCESSFULL COUNT")
return
figure(11)
setFigure(nofilterMatrix/sum(nofilter(fields)),"NO FILTER")
figure(12)
setFigure(filter1Matrix/sum(filter1(fields)),"FILTER 1")
figure(13)
setFigure(filter3Matrix/sum(filter3(fields)),"FILTER 3")
figure(14)
setFigure(resol1Matrix/sum(resol1(fields)),"RESOLUTION 1Hz")
figure(15)
setFigure(resol5Matrix/sum(resol5(fields)),"RESOLUTION 5Hz")
figure(16)
setFigure(resol10Matrix/sum(resol10(fields)),"RESOLUTION 10Hz")
figure(17)
setFigure(resol20Matrix/sum(resol20(fields)),"RESOLUTION 20Hz")
figure(18)
setFigure(resol50Matrix/sum(resol50(fields)),"RESOLUTION 50Hz")
figure(19)
setFigure(segm1Matrix/sum(segm1(fields)),"SEGMENT SIZE 1")
figure(20)
setFigure(segm2Matrix/sum(segm2(fields)),"SEGMENT SIZE 2")
figure(21)
setFigure(segm5Matrix/sum(segm5(fields)),"SEGMENT SIZE 5")
figure(22)
setFigure(segm10Matrix/sum(segm10(fields)),"SEGMENT SIZE 10")
figure(23)
setFigure(coher1Matrix/sum(coher1(fields)),"COHERENT FRACTIONS 1")
figure(24)
setFigure(coher3Matrix/sum(coher3(fields)),"COHERENT FRACTIONS 3")
figure(25)
setFigure(coher5Matrix/sum(coher5(fields)),"COHERENT FRACTIONS 5")


function setFigure(data,mtitle)
    load("reports\colorMAP1.mat")
    image(flipud(data),"CDataMapping","scaled")
    
    xt = [1:size(data,2)];
    yt = [1:size(data,1)];
    xtlbl = ["0" "2.5m" "7.5m" "17.5m" "37.5m" ...
             "77.5m" "157.5m" "317.5m" "637.5m" "637.5m" ...
             "1.2825" "1.6575" "2.0325" "3.155" "4.055"];
    %ytlbl = ["-14" "-18" "-22" "-26" "-30" ...
    %         "-32" "-34" "-36" "-38" "-40"...
    %         "-40.5" "-41" "-41.5" "-42" "-42.5"];
    ytlbl = ["46" "42" "38" "34" "30" ...
         "28" "26" "24" "22" "20"...
         "19.5" "19" "18.5" "18" "17.5"];
    set(gca, 'XTick',xt, 'XTickLabel',xtlbl)
    set(gca, 'YTick',yt, 'YTickLabel',ytlbl)
    xlabel("Doppler shift [Hz/sym]",Interpreter="latex")
    %ylabel("Signal/Noise [dB]",Interpreter="latex")
    ylabel("Received C/N0 [dB~Hz]",Interpreter="latex")
    title(mtitle)
    colormap(cc115); caxis([0 16]);
    colorbar('v',"XLim",[0 16],"XTick",[1 3 7 15],"XTickLabel",["RUNNED","ACQUIRED","TRACKED","DECODED"]);
end
function res = sumSub(final, new)
    newd = size(final)-size(new);
    newd = max(0, newd);
    padNew = padarray(new,newd,0,"post");
    padNew = padNew(1:size(final,1), 1:size(final,2));
    res = final + padNew;
end
function [res, ids, count] = maxSub(final, ids, count, new, newId)
    newd = size(final)-size(new);
    newd = max(0, newd);
    padNew = padarray(new,newd,0,"post");
    padNew = padNew(1:size(final,1), 1:size(final,2));
    for ab=1:numel(padNew)
        if padNew(ab) == 15
            count(ab) = count(ab) + 1;
            ids(ab) = newId;
        end
    end
    res = max(final, padNew);
end