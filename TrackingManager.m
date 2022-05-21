%
% First implementation: Lorenzo Borsoi & Davide Tomasella
% Review and Testing: Davide Tomasella & Lorenzo Borsoi
%
classdef TrackingManager < handle
    %Demodulator handles...

    properties 
        fsampling
        PRNsequence
        nPRN_x_Symbol
        %tracking evolution
        evolution struct
        currentStep
    end

    methods
        function obj = TrackingManager()
            %TrackingManager constructor
            %   Create TrackingManager class
            
            obj.evolution = obj.getNewEvolutionStepStruct();
            obj.currentStep = 0;
        end

        function obj = configCorrelatorValues(obj,fSampling, nPRN_x_Symbol, PRNcode)
            obj.fsampling = fSampling;
            obj.nPRN_x_Symbol = nPRN_x_Symbol;
            obj.PRNsequence = 2 * PRNcode - 1;   
        end

        function [decSymbols, idShift, idDoppler] = decodeOptimumShift(obj, inSamples, ...
                                                                       segmentSize, ...
                                                                       shifts_delayPRN, ...
                                                                       shifts_nSamples_x_symbolPeriod, ...
                                                                       shifts_nSamples_x_chipPeriod, ...
                                                                       coherenceFraction)    
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
            corrI = obj.normMultiply(shiftedPRNsampled, mySamples(1, :));
            corrQ = obj.normMultiply(shiftedPRNsampled, mySamples(2, :));
            
            %coherent integration, over the rows there are the coherent sums
            coherentCorrI = obj.sumOverCoherentFraction(corrI, shifts_delayPRN, ...
                                                        shifts_nSamples_x_symbolPeriod, ...
                                                        coherenceFraction, segmentSize);
            coherentCorrQ = obj.sumOverCoherentFraction(corrQ, shifts_delayPRN, ...
                                                        shifts_nSamples_x_symbolPeriod, ...
                                                        coherenceFraction, segmentSize);

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

            %complete correlation over symbols
            %TODO channel inversion? linear estimator?
            %NOTE: here we have to compensate the envelope rotation from 
            %      the start of the scenario (estimation during the 
            %      acquisition) before decoding the symbols
            bestCorrI = obj.sumFractionsOverSymbols(coherentCorrI(idMax, :), coherenceFraction);
            bestCorrQ = obj.sumFractionsOverSymbols(coherentCorrQ(idMax, :), coherenceFraction);
            
            %decoding
            decSymbols = (2 * (bestCorrI > 0) - 1)'; % column vector of decoded symbols +1,-1            
            
            figure(301)
            if obj.currentStep > 1
                hold on
            end
            set(gca,"ColorScale",'linear')
            %N.B.: X=delay,Y=doppler
            shading interp
            surf(obj.evolution(obj.currentStep).axis_delayPRN, ...
                 obj.evolution(obj.currentStep).axis_chipPeriod, ...
                 obj.evolution(obj.currentStep).trackingPeak, ...
                 'EdgeColor', 'none','FaceAlpha',0.4)
            hold off
            pause(0.3)
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

        function corr = normMultiply(~, shiftedPRNsampled, mySamples)
            corr = shiftedPRNsampled .* mySamples; % ./ sum(abs(shiftedPRNsampled) > 0, 2);
        end

        function coherentCorr = sumOverCoherentFraction(~, corr, shifts_delayPRN, ...
                                                        shifts_nSamples_x_symbolPeriod, ...
                                                        coherenceFraction, windowSize)
            nCoherentSegments = windowSize / coherenceFraction;
            maxOffset = max(shifts_delayPRN);
            optionLength = size(corr, 1);
            corrLength = size(corr, 2) - maxOffset;
            %row vector of coherent sum, option in each row
            coherentCorr = zeros(optionLength, nCoherentSegments);
            for cLine = 0:optionLength - 1 %cycle on columns
                symPeriod = shifts_nSamples_x_symbolPeriod(1 + mod(cLine, size(shifts_nSamples_x_symbolPeriod, 1)));
                symOffset = shifts_delayPRN(1 + fix(cLine / size(shifts_nSamples_x_symbolPeriod, 1)));
                %DT$ TODO check this approximation
                cohSamp = int32(symPeriod * coherenceFraction);
                cLine_padded = zeros(int32(corrLength / cohSamp + 0.5) * cohSamp,1);
                cLine_padded(1:corrLength) = corr(1 + cLine, 1 + symOffset:end - (maxOffset - symOffset))';
                corrSum = sum(reshape(cLine_padded, cohSamp, []), 1);
                coherentCorr(1 + cLine,:) = corrSum(1:nCoherentSegments);
            end
        end

        function bestCorr = sumFractionsOverSymbols(~, bestCoherentCorr, coherenceFraction)
            %columns of coherent fractions for each symbol
            tmpCorr = reshape(bestCoherentCorr, 1 / coherenceFraction, []);
            bestCorr = sum(tmpCorr, 1);
        end

        function stru = getNewEvolutionStepStruct(~)
            stru = struct("axis_delayPRN",[],"axis_chipPeriod",[], ...
                          "trackingPeak",[],"idDoppler",0,"idShift",0);                      
        end
     
    end
end