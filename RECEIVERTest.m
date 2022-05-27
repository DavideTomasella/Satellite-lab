classdef RECEIVERTest < matlab.unittest.TestCase
    
    properties (TestParameter)
        addDopplerID = {1,2,3,4,5,6,7,8,9,10}
        powerReduceID = {10,9,8,7,6,5,4,3,2,1}
        filterBandMultiplier = {0} %, 1 ,3
        % [1e2, 1e3, 5e3] / 100 = [1, 10, 50]
        reducedMaxDoppler = {1e2} %, 1e3, 5e3
        ppSegmentSize = {1, 2, 5, 10}
        nCoherentFractions = {1, 3, 5}
    end

    properties
        myREPORT = "REPORT_DT.mat"
        % NEW TEST FILE
        %TEST = "";
        %save("REPORT_DT.mat", "TEST");
        addDoppler = [    0,   0.2,   0.6,   1.4,     3, ...
                        6.2,  12.6,  25.4,    51, 102.6];
        powerReduce = [6.64,  6.31,  5.98,  5.65,  5.32, ...
                       4.98,  4.32,  3.65,  2.99,  2.32]
    end
    methods(TestClassSetup)
        % Shared setup for the entire test class
    end

    methods(TestMethodSetup)
        % Setup for each test
    end

    methods(Test)
        % Test methods

        function globalTest(testCase, ...
                            filterBandMultiplier, reducedMaxDoppler, ...
                            ppSegmentSize, nCoherentFractions, ...
                            addDopplerID, powerReduceID)
            mPowerReduce = testCase.powerReduce(powerReduceID);
            mAddDoppler = testCase.addDoppler(addDopplerID);

            try 
                REPORT = load(testCase.myREPORT);
                tmp = REPORT.(strcat("f"+num2str(filterBandMultiplier),"d"+num2str(reducedMaxDoppler), ...
                                "s"+num2str(ppSegmentSize),"c"+num2str(nCoherentFractions)));
                if tmp(powerReduceID,addDopplerID) > 0
                    %testCase.assumeTrue(true, "Test already run")
                    disp("Skipped test")
                    return
                end
                clear REPORT
                clear tmp
            catch
                if exist("REPORT","var")
                    clear REPORT
                end
                if exist("tmp","var")
                    clear tmp
                end
            end

            %% SIGNAL GENERATION
            tic
            %clearvars
            DEBUG = false;
            %powerReduce = 6
            %addDoppler = 0.3
            dopplerStart = 45.52;
            dopplerEnd = dopplerStart + mAddDoppler;
            signalDirectory = "binData/tempTest";
            testSignalName =  "T_tempTest";
            TESTPARS = struct("name",testSignalName, ...
                              "dstart",dopplerStart, ...
                              "dend",dopplerEnd, ...
                              "deveryChip",true, ...
                              "dmode",1, ...
                              "envelope",1, ...
                              "envelopePhase",0,...
                              "inNoise",0, ...
                              "outNoise",1, ...
                              "outNoise_Length",1223, ...
                              "outBits_Length",7, ...
                              "outPRN_Length",45, ...
                              "powerReduce",mPowerReduce);
            generateTestSignal
            toc
            pause(0.2)
            
            %%TEST RUN
            %clearvars
            tic
            DEBUG = false;
            newResult = false;
            outDirectory = "outData/tempTest";
            settingsName = "in0.json";
            resultsName = strcat("MM_test","_",num2str(mPowerReduce),"_",num2str(mAddDoppler), ...
                                 "_f",num2str(filterBandMultiplier),"_d",num2str(reducedMaxDoppler), ...
                                 "_s",num2str(ppSegmentSize),"_c",num2str(nCoherentFractions),".json");
            signalDirectory = "binData/tempTest";
            testSignalName =  "T_tempTest.bin";
            
            %filterBandMultiplier = 1
            %reduceMaxDoppler = 100
            nTestedDoppler = 100;
            centralDoppler = 50;
            thresholdSTD = 1; %correct threshold for tracking = 4
            %ppSegmentSize = 10
            %nCoherentFractions = 1
            
            RECEIVER
            toc
            status = uint8(1 + ...
                     2^1 * inout.results.ACQUISITION_OK + ...
                     2^2 * inout.results.TRACKING_OK + ...
                     2^3 * inout.results.DEMODULATION_OK);
            REPORT = load(testCase.myREPORT);
            REPORT = setfield(REPORT,"TREE","f"+num2str(filterBandMultiplier),"d"+num2str(reducedMaxDoppler), ...
                                "s"+num2str(ppSegmentSize),"c"+num2str(nCoherentFractions), ...
                                {powerReduceID,addDopplerID}, ...
                                status);
            REPORT = setfield(REPORT,strcat("f"+num2str(filterBandMultiplier),"d"+num2str(reducedMaxDoppler), ...
                                "s"+num2str(ppSegmentSize),"c"+num2str(nCoherentFractions)), ...
                                {powerReduceID,addDopplerID}, ...
                                status);
            save(testCase.myREPORT,'-struct',"REPORT",'-append');
            %TEST = "";
            %save("REPORT.mat","TEST");
            % MERGE TEST FILES
            %REPORT=load('REPORT_DT.mat')
            %save("REPORT.mat",'-struct',"REPORT",'-append');
        end
    end

end