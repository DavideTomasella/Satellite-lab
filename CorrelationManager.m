%
% First implementation: Gabriele Andreetta & Davide Tomasella (Acquisition algorithm)
%                       Davide Tomasella (Correlation peak update/tracking)
% Review and Testing:
%
classdef CorrelationManager < handle
    %Correlator handles...
    
    properties (SetAccess=private, GetAccess=public)
        fSampling
        nPRN_x_Symbol
        nChip_x_PRN
        PRNsequence
        SYNCsequence
        txChipRate
        %txSymbolRate @NOT_USED
        %search dimensions
        m_dopplerFreqs
        m_timeDelays;
        %corr results
        searchResults
        axis_delay
        axis_doppler
    end
    properties (SetAccess=public, GetAccess=public)
        %output real properties
        fDoppler {single}
        startingSample {uint32}
        %output virtual properties
        symbolPeriod {single}
        chipPeriod {single}
        nSamples_x_symbolPeriod {single}
        nSamples_x_chipPeriod {single}
        startingTime {single}
        initialPhase {single}
    end
        
    methods
        function obj = CorrelationManager()
            %CorrelationManager empty, no definitions required
            obj.searchResults = obj.getDefaultSearchResults();
        end        
        
        function obj = configCorrelatorMatrix(obj, fSampling, nPRN_x_Symbol, ...
                                              nChip_x_PRN, PRNcode, chipRate, ...
                                              syncPattern)
            %configCorrelatorMatrix configure parameters for calculation of
            %correlation matrix for signal acquisition
            %   fSampling: sampling frequencies
            %   nPRN_x_Symbol: n. of repeated PRN for symbol period
            %   nChip_x_PRN: length of PRN sequence
            %   PRNcode: binary vector with PRN sequence
            %   chipRate: transmitter chip rate (without doppler)
            %   syncPattern: string with synchronization pattern

            obj.fSampling = fSampling;
            obj.nPRN_x_Symbol = nPRN_x_Symbol;
            obj.PRNsequence = 2 * PRNcode - 1;
            if length(PRNcode) ~= nChip_x_PRN
                disp("Error, PRN code is incompatible with modulation parameters")
            end
            obj.nChip_x_PRN = nChip_x_PRN;
            obj.txChipRate = chipRate;
            %obj.txSymbolRate = chipRate / nChip_x_PRN / nPRN_x_Symbol; @NOT_USED
            syncPattern = decimalToBinaryVector(bin2dec(syncPattern), ...
                                                length(char(syncPattern)));
            obj.SYNCsequence = 2 * syncPattern - 1;
        end        
        
        function [maxMatrix, meanMatrix, squareMatrix] = ...
                                    calcCorrelationMatrix(obj, IQsamples, dimMatrix, ...
                                                            maxDoppler, fresolution, currentSample)
            %calcCorrelationMatrix handles the creation of the correlation
            %matrix (doppler frequencies and time delays)
            % ##Code adapted from:
            % ##source: http://maxwell.sze.hu/~ungert/Radiorendszerek_satlab/Segedanyagok/Ajanlott_irodalom/matlab_simulations_for_radar_systems_design_(cuppy).pdf
            % ##page: 247/686
            %   IQsamples: column vector with input samples
            %   dimMatrix: dimension NxN of output correlation matrix
            %   maxDoppler: maximum doppler
            %   fresolution: doppler frequency step size
            %   currentSample: delay (in samples) of the considered signal

            delay_redFactor = fix(length(IQsamples) / 2 / dimMatrix);
            freq_redFactor = ceil((2 * maxDoppler / fresolution) / dimMatrix);
            delay_redFactor = max(delay_redFactor, 1);
            freq_redFactor = max(freq_redFactor, 1);
            Ndelays = delay_redFactor * dimMatrix;
            Nfrequencies = freq_redFactor * dimMatrix;
            Nsamples = Ndelays * 2;
            
            obj.m_timeDelays = (0:Nsamples-1)' / obj.fSampling;
            obj.m_dopplerFreqs = linspace(-maxDoppler, maxDoppler, Nfrequencies);
            %matrices for plot results
            maxMatrix = zeros(dimMatrix);
            meanMatrix = zeros(dimMatrix);
            squareMatrix = zeros(dimMatrix);
            obj.axis_delay = obj.m_timeDelays(1:delay_redFactor:Ndelays);
            obj.axis_doppler = obj.m_dopplerFreqs(1:freq_redFactor:Nfrequencies);
            
            %NOTE: if false find correlation only for a symbol (some PRNs)!
            useSyncPattern = true;

            %reset searching results
            obj.searchResults = obj.getDefaultSearchResults();
            
            %compute input fft
            input_fft = fft(IQsamples(:, 1) + 1i * IQsamples(:, 2), Nsamples);
            
            %cycle over frequencies
            for h = 1:Nfrequencies

                %create interp PRN
                nChipPeriod = obj.fSampling / (obj.txChipRate + obj.m_dopplerFreqs(h));
                PRNsampled = obj.getPRNFromChipPeriod(nChipPeriod, Nsamples, useSyncPattern);

                % taken from ambiguity function
                PRNsampled_fft = fft(PRNsampled.* ...
                                     exp(1i * 2 * pi * obj.m_dopplerFreqs(h) * ...
                                         (obj.m_timeDelays + currentSample / obj.fSampling)), ...
                                     Nsamples);
 
                %GA$ we need the phase for the downconverted and the
                %demodulation of the symbols (phase envelope tracking)
                %delayCorrelation = abs(ifft(PRNsampled_fft .* conj(input_fft))) / Nsamples;
                complexCorrelation = ifft(PRNsampled_fft .* conj(input_fft)) / Nsamples;
                delayCorrelation = abs(complexCorrelation);
                
                % calculate and save peak precise position and its phase
                [max_delayCorrelation, maxIndex] = max(delayCorrelation, [], 1);
                phaseCorrelation = angle(complexCorrelation(maxIndex));
                if max_delayCorrelation > obj.searchResults.maxPeak
                    obj.searchResults.maxPeak = max_delayCorrelation;
                    obj.searchResults.idStartTime = maxIndex;
                    obj.searchResults.idDopplerShift = h;
                    obj.searchResults.phase = phaseCorrelation;
                end
                %GA$ complex mean and mean squared values of correlation
                obj.searchResults.mean = obj.searchResults.mean + sum(complexCorrelation) / Nsamples / Nsamples;
                obj.searchResults.meanSquare = obj.searchResults.meanSquare +...
                                            sum(complexCorrelation.^2) / Nsamples / Nsamples;
                
                %save reduced matrix, lineCorrelation is a column
                reducedLine = max(reshape(delayCorrelation, delay_redFactor, []), [], 1); %row vector
                maxMatrix(:, ceil(h / freq_redFactor)) = max(maxMatrix(:, ceil(h / freq_redFactor)), ...
                                                           reducedLine(1:dimMatrix)');
                reducedMeanLine = sum(reshape(delayCorrelation, delay_redFactor, []), 1)...
                                    / delay_redFactor / freq_redFactor;
                meanMatrix(:, ceil(h / freq_redFactor)) =...
                            sum([meanMatrix(:, ceil(h / freq_redFactor)) reducedMeanLine(1:dimMatrix)'],2);
                reducedSquaredLine = sum(reshape(delayCorrelation, delay_redFactor, []).^2, 1)...
                                    / delay_redFactor / freq_redFactor;
                squareMatrix(:, ceil(h / freq_redFactor)) =...
                            sum([squareMatrix(:, ceil(h / freq_redFactor)) reducedSquaredLine(1:dimMatrix)'],2);
                if mod(h, 10 * freq_redFactor) == 0
                    sprintf("Completed %0.1f%%", h / Nfrequencies * 100)
                    figure(201)
                    set(gca, "ColorScale", 'log')
                    image(obj.axis_doppler, obj.axis_delay, maxMatrix, 'CDataMapping', 'scaled')
                    pause(1)
                end
            end
            %GA$ absolute value of mean, mean squared, mean matrix and squared
            % matrix (otherwise are complex values)
            obj.searchResults.mean = abs(obj.searchResults.mean);
            obj.searchResults.meanSquare = abs(obj.searchResults.meanSquare);
            meanMatrix = abs(meanMatrix);
            squareMatrix = abs(squareMatrix);
            %obj.searchResults
            
            %close(201)
            figure(202)
            set(gca,"ColorScale",'linear')
            surf(obj.axis_doppler, obj.axis_delay, maxMatrix, 'EdgeColor', 'none')
            pause(1)
        end

        function PRNsampled = getPRNFromChipPeriod(obj, nSamples_x_chipPeriod, ...
                                                   Nsamples, useSyncPattern)
            %maxLength = nSamples_x_chipPeriod * length(obj.PRNsequence) * ...
            %            obj.nPRN_x_Symbol;
            %if useSyncPattern
            %    maxLength = maxLength * length(obj.syncPattern);
            %end
            PRNsampled = zeros(Nsamples, 1); %ceil approx

            PRNsymbol = repmat(obj.PRNsequence, 1, obj.nPRN_x_Symbol);
            if useSyncPattern
                PRNsync = reshape(obj.SYNCsequence .* PRNsymbol', [], 1);
            else
                PRNsync = PRNsymbol';  
            end
            newTime = 0:1 / nSamples_x_chipPeriod:length(PRNsync);
            PRNsampled(1:length(newTime)) = interp1((0:length(PRNsync))', [PRNsync; PRNsync(1)], ...
                                                    newTime', "previous"); %upsampling
        end

        function results = getDefaultSearchResults(~)
            results = struct("maxPeak", 0, ... %maximum is abolute value
                             "idStartTime", 0, ...
                             "idDopplerShift", 0, ...
                             "mean", 0, ...
                             "meanSquare", 0, ...
                             "phase", 0); 
        end
        
        %                           %
        %     DAVIDE TOMASELLA      %
        %  Functions for tracking.  %
        %  The parameter updates    %
        %  after each step.         %
        %                           %

        %DT$ already corrected
        function isAcquired = ifAcquiredFindCorrelationPeak(obj, thresholdSTD, ...
                                                            outInterface)
            %ifAcquiredFindCorrelationPeak handles the decision if a peak
            %is found and the saving of the correlation peak parameters
            %The correlation acquisition results are contained in obj.searchResults
            %   thresholdSTD: n. STD over average to identify a correlation peak
            %   outInterface: outInterface.results contains the results

            % define dynamic threshold
            std1 = sqrt(obj.searchResults.meanSquare - obj.searchResults.mean ^ 2);
            thresh = obj.searchResults.mean + thresholdSTD * std1;
            if (obj.searchResults.maxPeak > 0 && obj.searchResults.maxPeak > thresh)
                %get bet doppler and delay
                obj.fDoppler = obj.m_dopplerFreqs(obj.searchResults.idDopplerShift);
                obj.startingTime = obj.m_timeDelays(obj.searchResults.idStartTime);
                obj.initialPhase = obj.searchResults.phase;
                %output saving
                outInterface.results.estimatedDopplerStart = obj.fDoppler;
                outInterface.results.estimatedDelay = obj.startingTime;
                outInterface.results.estimatedPhase = obj.initialPhase;
                %return successfull acquisition flag
                isAcquired = true;
            else
                %sync pattern not found not found
                isAcquired = false;
            end
        end
        
        %DT$ already corrected
        function obj = updateCorrelationPeak(obj, new_SamplesChipPeriod, advancement_startingSample)
            %updateCorrelationPeak handles the update of the peak position
            %during the tracking phase and keep updated the parameters
            %   new_SamplesSymbolPeriod: n. samples in new symbol period
            %   advancement_startingSample: n. samples to advance
            if new_SamplesChipPeriod<=0
                warning("Error. The number of samples per symbol period cannot be a negative value.")
            end
            if advancement_startingSample<=0
                warning("Error. The advancement size cannot be a negative value.")
            end
            %DT$ not needed thanks to dynamic properties
            %newChipPeriod = new_symbolPeriod / obj.nPRN_x_Symbol / obj.nChip_x_PRN;
            %obj.fDoppler = (1 / newChipPeriod) - obj.fModulation;
            %oldDoppler = obj.fDoppler;
            newChipPeriod = new_SamplesChipPeriod / obj.fSampling;
            advanceTime = single(advancement_startingSample) / obj.fSampling;
            %advancePhase = 2 * pi / obj.chipPeriod * obj.startingTime - ...
            %               2 * pi / newChipPeriod * (obj.startingTime + advanceTime);
            advancePhase = 2 * pi * obj.fDoppler * advanceTime; % + advancement_phase
            obj.chipPeriod = newChipPeriod;
            obj.startingSample = obj.startingSample + uint32(advancement_startingSample);
            obj.initialPhase = obj.initialPhase + advancePhase;
        end

        %%%%%%%%% REAL PROPERTIES %%%%%%%%%

        %function set.fDoppler

        %function get.fDoppler 

        function set.startingSample(obj, iStartingSample)
            obj.startingSample = uint32(iStartingSample);
        end
       
        %function get.startingSample

        function set.initialPhase(obj, iInitialPhase)
            obj.initialPhase = mod(iInitialPhase, 2 * pi);
        end

        %function get.initialPhase

        %%%%%%%% VIRUAL PROPERTIES %%%%%%%%

        function oSamplesSymbolPeriod = get.nSamples_x_symbolPeriod(obj)
            oSamplesSymbolPeriod = obj.symbolPeriod * obj.fSampling;
        end

        function set.nSamples_x_symbolPeriod(obj, iSamplesSymbolPeriod)
            obj.symbolPeriod = iSamplesSymbolPeriod / obj.fSampling;
        end

        function oSymbolPeriod = get.symbolPeriod(obj)
            oSymbolPeriod = obj.chipPeriod * obj.nChip_x_PRN * obj.nPRN_x_Symbol;
        end

        function set.symbolPeriod(obj, iSymbolPeriod)
            obj.chipPeriod = iSymbolPeriod / obj.nPRN_x_Symbol / obj.nChip_x_PRN;
        end

        function oSamplesChipPeriod = get.nSamples_x_chipPeriod(obj)
            oSamplesChipPeriod = obj.chipPeriod * obj.fSampling;
        end

        function set.nSamples_x_chipPeriod(obj, iSamplesChipPeriod)
            obj.chipPeriod = iSamplesChipPeriod / obj.fSampling;
        end

        function oChipPeriod = get.chipPeriod(obj)
            oChipPeriod = 1 / (obj.txChipRate + obj.fDoppler);
        end

        function set.chipPeriod(obj, iChipPeriod)
            obj.fDoppler = (1 / iChipPeriod) - obj.txChipRate;
        end

        function oStartingTime = get.startingTime(obj)
            oStartingTime = single(obj.startingSample) / obj.fSampling;
        end

        function set.startingTime(obj, iStartingTime)
            obj.startingSample = uint32(iStartingTime * obj.fSampling);
        end

    end
end
