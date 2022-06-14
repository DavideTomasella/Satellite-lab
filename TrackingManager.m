%
% First implementation: Lorenzo Borsoi & Davide Tomasella
% Review and Testing: Davide Tomasella & Lorenzo Borsoi
%
classdef TrackingManager < handle
    %Demodulator handles...

    properties (SetAccess=private, GetAccess=public)
        fsampling
        PRNsequence
        nPRN_x_Symbol
        maxSignalValue
        %PLL-MSEtracking
        %SYNCsequence @NOT_USED
        MMSEphase
        MMSEalpha
        %tracking evolution
        evolution struct
        currentStep
        plotID1
        plotID2
    end

    properties (SetAccess=public, GetAccess=public)
        DEBUG
    end

    methods
        function obj = TrackingManager(DEBUG)
            %TrackingManager constructor
            %   Create TrackingManager class
            if nargin < 1
                DEBUG = false;
            end

            obj.DEBUG = DEBUG;            
            obj.evolution = obj.getNewEvolutionStepStruct();
            obj.currentStep = 0;
        end

        function obj = configCorrelatorValues(obj,fSampling, nPRN_x_Symbol,  ...
                                              PRNcode, MMSEalpha, quantizationBits)
            obj.fsampling = fSampling;
            obj.nPRN_x_Symbol = nPRN_x_Symbol;
            obj.PRNsequence = 2 * PRNcode - 1;   
            %syncPattern = decimalToBinaryVector(bin2dec(syncPattern), ...
            %                                    length(char(syncPattern)));
            %obj.SYNCsequence = 2 * syncPattern - 1;
            %tracking through atan: insensitive to pi shift
            obj.MMSEphase = 0;
            obj.MMSEalpha = MMSEalpha;
            obj.maxSignalValue = 2 ^ (quantizationBits - 1);
        end

        function [decSymbols, idShift, idDoppler] = decodeOptimumShift(obj, inSamples, ...
                                                                       segmentSize, ...
                                                                       shifts_delayPRN, ...
                                                                       shifts_nSamples_x_symbolPeriod, ...
                                                                       shifts_nSamples_x_chipPeriod, ...
                                                                       nCoherentFractions, ...
                                                                       outInterface)    
            %windowSize: numero di simboli
            %nchips: numero di chips che compongono la PNR sequence
            %numero di PNR sequence in un simbolo: symbolPeriod/(chipPeriod*nchips)
            %numero di PNR sequence in M messaggi: windowSize*symbolPeriod/(chipPeriod*nchips)
            %PNRsequence va upsampled con fsampling in PNRcodes 1xN

            %Reshape input to have predetermined vector
            if size(inSamples, 1) < size(inSamples, 2)
                inSamples = inSamples'; %column vector
            end
            if size(shifts_nSamples_x_chipPeriod, 1) < size(shifts_nSamples_x_chipPeriod, 2)
                shifts_nSamples_x_chipPeriod = shifts_nSamples_x_chipPeriod'; %column vector
            end
            if size(shifts_nSamples_x_symbolPeriod, 1) < size(shifts_nSamples_x_symbolPeriod, 2)
                shifts_nSamples_x_symbolPeriod = shifts_nSamples_x_symbolPeriod'; %column vector
            end
            if size(shifts_delayPRN, 2) < size(shifts_delayPRN, 1)
                shifts_delayPRN = shifts_delayPRN'; %row vector
            end
            nCoherentFractions = ceil(nCoherentFractions); %must be integer to decode symbols
            
            %Initialization evolution monitoting
            %NOTE: (DT) shifts_delayPRN changes, I want to save the input
            %           one instead of converting it explicitly -(end+1)/2
            obj.currentStep = obj.currentStep + 1;
            obj.evolution(obj.currentStep) = obj.getNewEvolutionStepStruct();
            obj.evolution(obj.currentStep).axis_chipPeriod = shifts_nSamples_x_chipPeriod;
            obj.evolution(obj.currentStep).axis_delayPRN = shifts_delayPRN;

            %create PRN with different length due to different dopplers
            PRNsampled = obj.getPRNFromChipPeriods(shifts_nSamples_x_chipPeriod, segmentSize);

            %adapt the length of the samples to the PRN, horizontal vector
            [mySamples, shifts_delayPRN] = obj.adaptSamplesToPRNlength(inSamples, shifts_delayPRN, ...
                                                                       size(PRNsampled, 2));
            
            %add different delays to the PRN, sampled PRN with shifts in time and frequency
            shiftedPRNsampled = obj.createShiftedPRN(PRNsampled, shifts_delayPRN);
            %plot(shiftedPRNsampled(1:2,end-1e4:1:end)')
            
            %in-phase & quadrature multicorrelation
            corrI = obj.normMultiply(shiftedPRNsampled, mySamples(1, :), segmentSize);
            corrQ = obj.normMultiply(shiftedPRNsampled, mySamples(2, :), segmentSize);
            
            %coherent integration, over the rows there are the coherent sums
            coherentCorrI = obj.sumOverCoherentFraction(corrI, shifts_delayPRN, ...
                                                        shifts_nSamples_x_symbolPeriod, ...
                                                        nCoherentFractions, segmentSize);
            coherentCorrQ = obj.sumOverCoherentFraction(corrQ, shifts_delayPRN, ...
                                                        shifts_nSamples_x_symbolPeriod, ...
                                                        nCoherentFractions, segmentSize);

            %non-coherent integration, column vector
            noncoherentCorr = sum(coherentCorrI .^ 2 + coherentCorrQ .^ 2, 2);
            obj.evolution(obj.currentStep).trackingPeak = reshape(noncoherentCorr, ...
                                                                  size(shifts_nSamples_x_chipPeriod, 1), ...
                                                                  size(shifts_delayPRN, 2));

            %find max
            [~, idMax] = max(noncoherentCorr,[], 1);
            [idDoppler, idShift] = ind2sub([size(shifts_nSamples_x_chipPeriod, 1) size(shifts_delayPRN, 2)], idMax);         
            %NOTE: DOT product: sommatoria[(Ie-Il)*Ip] - sommatoria[(Qe-Ql)*qp],
            %if <0 detector is late, if >0 too early
            obj.evolution(obj.currentStep).idDoppler = idDoppler;
            obj.evolution(obj.currentStep).idShift = idShift;
            
            %if obj.DEBUG
            %    figure
            %    j1=sum(coherentCorrI,2)
            %    j2=sum(coherentCorrQ,2)
            %    plot3(1:length(j1),j1,j2,"o-")
            %    grid on
            %end
            
            %compensate error in downconvertion
            rotatedCorr = obj.compensateResidualEnvelope(coherentCorrI(idMax, :), coherentCorrQ(idMax, :));

            %coherent sum over symbols
            bestCorr = obj.sumFractionsOverSymbols(rotatedCorr, nCoherentFractions);

            %decoding
            decSymbols = (2 * (real(bestCorr) > 0) - 1)'; % column vector of decoded symbols +1,-1            
            trackOK = sum(abs(real(bestCorr) ./ imag(bestCorr))) > segmentSize * nCoherentFractions; 
            outInterface.results.TRACKING_OK = outInterface.results.TRACKING_OK & trackOK;

            sprintf("Tracker step %d: demodulated %d symbols.", obj.currentStep, obj.currentStep * segmentSize)
            if obj.DEBUG
                disp(sum(abs(real(bestCorr) ./ imag(bestCorr))) / segmentSize / nCoherentFractions)
                obj.plotTrackingPeak();
                obj.plotDemodulatedSymbols(bestCorr, segmentSize);   
                pause(0.3)             
            end

        end

        %                        %
        % PRIVATE METHODS        %
        %                        %

        function PRNsampled = getPRNFromChipPeriods(obj, shifts_nSamples_x_chipPeriod, windowSize)
            optionDopplerLength = size(shifts_nSamples_x_chipPeriod, 1);
            maxLength = max(shifts_nSamples_x_chipPeriod * length(obj.PRNsequence) * ...
                            obj.nPRN_x_Symbol * windowSize);
            PRNsampled = zeros(optionDopplerLength, int32(maxLength + 0.5)); %ceil approx

            for period = 1:optionDopplerLength
                repSequence = repmat(obj.PRNsequence, 1, windowSize * obj.nPRN_x_Symbol);                
                PRNinterp = interp1(0:length(repSequence), [repSequence repSequence(1)], ...
                    0:1 / shifts_nSamples_x_chipPeriod(period):length(repSequence), "previous"); %upsampling   
                repSize = size(PRNinterp,2);
                PRNsampled(period, 1:repSize) = PRNinterp;
            end
        end

        function  [mySamples, shifts_delayPRN] = adaptSamplesToPRNlength(~, filteredSamples, ...
                                                                         shifts_delayPRN, PRNlength)
            offset = -min(shifts_delayPRN);
            newLength = PRNlength - offset;
            if size(filteredSamples, 1) < newLength
                mySamples = [zeros(offset, 2); filteredSamples; zeros(newLength - size(filteredSamples, 1), ...
                                                    newLength)]';
            else
                mySamples = [zeros(offset, 2); filteredSamples(1:newLength, :)]';
            end
            shifts_delayPRN = shifts_delayPRN + offset;
        end

        function shiftedPRNsampled = createShiftedPRN(~, PRNsampled, shifts_delayPRN)
            optionDopplerLength = size(PRNsampled, 1);
            optionDelayLength = size(shifts_delayPRN, 2);
            shiftedPRNsampled = zeros(optionDelayLength * optionDopplerLength, ...
                                      size(PRNsampled, 2));
            for shiftID = 0:optionDelayLength - 1
                offset = shiftID * optionDopplerLength;
                %circ shift of PRN sequence
                shiftedPRNsampled(1 + offset:offset + optionDopplerLength, :) = ...
                            circshift(PRNsampled, shifts_delayPRN(1 + shiftID), 2);
            end
        end

        function corr = normMultiply(obj, shiftedPRNsampled, mySamples, segmentSize)
            corr = shiftedPRNsampled .* mySamples * ...
                            segmentSize / length(mySamples) / obj.maxSignalValue; 
            %normalization to 1 * nSymbols 
        end

        function coherentCorr = sumOverCoherentFraction(~, corr, shifts_delayPRN, ...
                                                        shifts_nSamples_x_symbolPeriod, ...
                                                        nCoherentFractions, windowSize)
            nCoherentSegments = windowSize * nCoherentFractions;
            maxOffset = max(shifts_delayPRN);
            optionLength = size(corr, 1);
            corrLength = size(corr, 2) - maxOffset;
            %row vector of coherent sum, option in each row
            coherentCorr = zeros(optionLength, nCoherentSegments);
            for cLine = 0:optionLength - 1 %cycle on columns
                symPeriod = shifts_nSamples_x_symbolPeriod(1 + mod(cLine, size(shifts_nSamples_x_symbolPeriod, 1)));
                symOffset = shifts_delayPRN(1 + fix(cLine / size(shifts_nSamples_x_symbolPeriod, 1)));
                %DT$ TODO check this approximation
                cohSamp = int32(symPeriod / nCoherentFractions);
                cLine_padded = zeros(int32(corrLength / cohSamp + 0.5) * cohSamp,1);
                cLine_padded(1:corrLength) = corr(1 + cLine, 1 + symOffset:end - (maxOffset - symOffset))';
                corrSum = sum(reshape(cLine_padded, cohSamp, []), 1);
                coherentCorr(1 + cLine,:) = corrSum(1:nCoherentSegments);
            end
        end

        function bestCorr = sumFractionsOverSymbols(~, rotatedCoherentCorr, nCoherentFractions)
            %columns of coherent fractions for each symbol
            tmpCorr = reshape(rotatedCoherentCorr, nCoherentFractions, []);
            bestCorr = sum(tmpCorr, 1);
        end

        function rotatedCorr = compensateResidualEnvelope(obj, bestCoherentCorrI, bestCoherentCorrQ)
            %estimate angle shift (estimated doppler - real doppler)
            %through MMSE algorithm -> second order filter = PLL
            preCompensatedCorrA = (bestCoherentCorrI + 1i * bestCoherentCorrQ) .* ...
                exp(-1i * obj.MMSEphase);
            preAngleA = atan(imag(preCompensatedCorrA) ./ real(preCompensatedCorrA));
            coherentAngle = [0 preAngleA];
            %correct the pi transitions due to -pi/2 pi/2 interval for atan
            %NOTE: since atan has period pi we neglect the symbol phase (0-pi)
            piT1 = find(diff(coherentAngle) > pi - pi/2);
            for t = piT1
                coherentAngle(1+t) = coherentAngle(1+t) - pi;
                %disp("inv+")
            end
            piT2 = find(diff(coherentAngle) < -pi + pi/2);
            for t = piT2
                coherentAngle(1+t) = coherentAngle(1+t) + pi;
                %disp("inv-")
            end
            %if obj.DEBUG
            %    figure(30)
            %    hold on
            %    plot(coherentAngle)
            %    hold off
            %end

            filteredAngle = filter([obj.MMSEalpha, 1 - obj.MMSEalpha], 1, ...
                                   coherentAngle);
            obj.MMSEphase = mod(obj.MMSEphase + filteredAngle(end) + pi, 2 * pi) - pi;
            %invert rotation ans sum over symbols
            rotatedCorr = preCompensatedCorrA .* ...
                            exp(-1i * filteredAngle(2:end));
        end


        function rotatedCorr = compensateResidualEnvelope2(obj, bestCoherentCorrI, bestCoherentCorrQ)
            %estimate angle shift (estimated doppler - real doppler)
            %through MMSE algorithm -> second order filter = PLL
            preCompensatedCorrA = (bestCoherentCorrI + 1i * bestCoherentCorrQ) .* ...
                            exp(-1i * obj.MMSEphase);
            preCompensatedCorrB = (bestCoherentCorrI + 1i * bestCoherentCorrQ) .* ...
                exp(-1i * (obj.MMSEphase - pi/4));
            preCompensatedCorrC = (bestCoherentCorrI + 1i * bestCoherentCorrQ) .* ...
                exp(-1i * (obj.MMSEphase + pi/4));
            coherentAngle = atan(bestCoherentCorrQ ./ bestCoherentCorrI);
            preAngleA = atan(imag(preCompensatedCorrA) ./ real(preCompensatedCorrA));
            preAngleB = atan(imag(preCompensatedCorrB) ./ real(preCompensatedCorrB)) - pi/4;
            preAngleC = atan(imag(preCompensatedCorrC) ./ real(preCompensatedCorrC)) + pi/4;
            preAngleD = atan2(imag(preCompensatedCorrA), real(preCompensatedCorrA));
            %plot(preAngleA)
            %hold on
            %plot(preAngleB)
            %plot(preAngleC)
            %plot(preAngleB-preAngleA)
            %plot(preAngleC-preAngleA)
            %hold off
            
            %find transitions
            upTransitions = (preAngleB-preAngleA < -pi + pi/2);
            downTransitions = (preAngleC-preAngleA > pi - pi/2);
            %accept the ones where the new value is closer to zero
            postAngleA = [0 preAngleA];
            for id = 1:length(preAngleA)
                newAngleA = postAngleA;
                if upTransitions(id)
                    newAngleA(id+1) = newAngleA(id+1) - pi;
                elseif downTransitions(id)
                    newAngleA(id+1) = newAngleA(id+1) + pi;
                end
                oldDiff = diff(postAngleA);
                newDiff = diff(newAngleA);
                previousSum = cumsum(postAngleA);
                oldDist = postAngleA(id+1) - previousSum(id) / id;
                newDist = newAngleA(id+1) - previousSum(id) / id;
                if abs(newDiff(id)) < pi / 4 
                    disp("diff")
                    postAngleA = newAngleA;
                elseif abs(newDist) < abs(oldDist)
                    postAngleA = newAngleA;
                end
            end
                
            %diffUpTransitions = abs(diff(preAngleA - pi * upTransitions)) - abs(diff(preAngleA));
            %diffDownTransitions = abs(diff(preAngleA + pi * downTransitions)) - abs(diff(preAngleA));
            %acceptedUpTransitions = upTransitions .* ([0 diffUpTransitions] < 0);
            %acceptedDownTransitions = downTransitions .* ([0 diffDownTransitions] < 0);
            %postAngleA = preAngleA + pi * acceptedDownTransitions - pi * acceptedUpTransitions;
            if false
            coherentAngle = [0 preAngleA];
            else
            coherentAngle = [obj.MMSEphase coherentAngle];
            end
            %correct the pi transitions due to -pi/2 pi/2 interval for atan
            %NOTE: since atan has period pi we neglect the symbol phase (0-pi)
            piT1 = find(diff(coherentAngle) > pi - pi/2);
            for t = piT1
                disp("inv+")
                coherentAngle(1+t:end) = coherentAngle(1+t:end) - pi;
            end
            piT2 = find(diff(coherentAngle) < -pi + pi/2);
            for t = piT2
                disp("inv-")
                coherentAngle(1+t:end) = coherentAngle(1+t:end) + pi;
            end
            if obj.DEBUG
                h=figure(303);
                movegui(h,"northeast")
                movegui(h,"center")
                hold off
                plot(preAngleA)
                hold on
                plot(preAngleD)
                plot(abs(preAngleD)<pi/2)
                %plot(pi * upTransitions)
                %plot(abs(preAngleB) < abs(preAngleA))
                %plot(-pi * downTransitions)
                %plot(abs(preAngleC) - abs(preAngleA))
                plot(postAngleA,"--")
                plot(abs([0 preAngleD]-postAngleA)<pi/2,".--")
            end
            %if obj.DEBUG
            %    figure(30)
            %    hold on
            %    plot(coherentAngle)
            %    hold off
            %end
            
            if true
            filteredAngle = filter([obj.MMSEalpha, 1 - obj.MMSEalpha], 1, ...
                                   coherentAngle);
            elseif false
            filteredAngle = filter([obj.MMSEalpha, 1 - obj.MMSEalpha], 1, ...
                                   postAngleA);
            else
            filteredAngle = filter([1 - obj.MMSEalpha, obj.MMSEalpha], 1, ...
                                   [0 preAngleA]);
            end
            
            if false
            %obj.MMSEphase = mod(filteredAngle(end) + pi, 2 * pi) - pi;
            %obj.MMSEphase = mod(obj.MMSEphase + filteredAngle(end-1) + pi, 2 * pi) - pi;
            obj.MMSEphase = obj.MMSEphase + filteredAngle(end);
            %invert rotation ans sum over symbols
            rotatedCorr = preCompensatedCorrA .* ...
                            exp(-1i * filteredAngle(2:end));
            else
            filteredAngle = filter([obj.MMSEalpha, 1 - obj.MMSEalpha], 1, ...
                                   coherentAngle);
            obj.MMSEphase = mod(filteredAngle(end) + pi, 2 * pi) - pi;
            %obj.MMSEphase = mod(obj.MMSEphase + filteredAngle(end-1) + pi, 2 * pi) - pi;
            %obj.MMSEphase = obj.MMSEphase + filteredAngle(end);
            %invert rotation ans sum over symbols
            rotatedCorr = (bestCoherentCorrI + 1i * bestCoherentCorrQ) .* ...
                            exp(-1i * filteredAngle(2:end));
            end

        end


        function stru = getNewEvolutionStepStruct(~)
            stru = struct("axis_delayPRN",[],"axis_chipPeriod",[], ...
                          "trackingPeak",[],"idDoppler",0,"idShift",0);                      
        end

        function plotTrackingPeak(obj)
            fh301 = figure(301);
            movegui(fh301,"south")
            if obj.currentStep == 1
                obj.plotID1 = plot3(0,0,0,".-", MarkerSize=8, Color=[0.7 0 0 1]);
                obj.plotID1.XDataSource='xPKK';
                obj.plotID1.YDataSource='yPKK';
                obj.plotID1.ZDataSource='zPKK';
            end
            hold on
            set(gca,"ColorScale",'linear')
            %N.B.: X=delay,Y=doppler
            shading interp
            surf(obj.evolution(obj.currentStep).axis_delayPRN, ...
                 obj.evolution(obj.currentStep).axis_chipPeriod, ...
                 obj.evolution(obj.currentStep).trackingPeak, ...
                 'EdgeColor', 'none','FaceAlpha',0.4)
            xPKK = obj.plotID1.XData;
            yPKK = obj.plotID1.YData;
            zPKK = obj.plotID1.ZData;
            xPKK(obj.currentStep) = obj.evolution(obj.currentStep).axis_delayPRN(obj.evolution(obj.currentStep).idShift);
            yPKK(obj.currentStep) = obj.evolution(obj.currentStep).axis_chipPeriod(obj.evolution(obj.currentStep).idDoppler);
            zPKK(obj.currentStep) = obj.evolution(obj.currentStep).trackingPeak(obj.evolution(obj.currentStep).idDoppler,...
                                                                                obj.evolution(obj.currentStep).idShift);

            hold off
            xlabel("Delay shift [tens of chip period]")
            ylabel("Frequency shift [samples per chip period]")
            zlabel("Correlation peak (normalized per symbol)")
            grid on
            view([27.5 50])
            refreshdata(fh301, 'caller')
        end

        function plotDemodulatedSymbols(obj, bestCorr, segmentSize)
            fh302 = figure(302);
            movegui(fh302,"southeast")
            if obj.currentStep == 1
                obj.plotID2 = plot(0,0,".:", MarkerSize=10, Color=[0 0 0.8 1]);
                %xlim([-1.1 1.1])
                %ylim([-1.1 1.1])
                obj.plotID2.XDataSource='xSYM';
                obj.plotID2.YDataSource='ySYM';
            end
            xSYM = obj.plotID2.XData;
            ySYM = obj.plotID2.YData;
            xSYM(1 + (obj.currentStep - 1) * segmentSize:obj.currentStep * segmentSize) = ...
                    real(bestCorr); % / sqrt(peakValue)
            ySYM(1 + (obj.currentStep - 1) * segmentSize:obj.currentStep * segmentSize) = ...
                    imag(bestCorr); % / sqrt(peakValue)
            xlabel("In-phase correlation (normalized per symbol)")
            ylabel("Quadrature correlation (normalized per symbol)")
            xlim([-1.1 1.1] .* max(abs(xSYM)))
            ylim([-2 2] .* max(abs(ySYM)))
            grid on
            refreshdata(fh302, 'caller')
        end
     
    end
end