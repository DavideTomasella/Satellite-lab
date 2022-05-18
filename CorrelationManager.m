%
% First implementation: Davide Tomasella & ... (Acquisition algorithm)
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
        txChipRate
        %txSymbolRate @NOT_USED
        syncPattern
        %search dimensions
        m_chipPeriods
        m_timeDelays

        maxCorrelation
        dopplerFrequencies
        delays;
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
    end
        
    methods
        function obj = CorrelationManager()
            %CorrelationManager empty, no definitions required
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
            obj.syncPattern = decimalToBinaryVector(bin2dec(syncPattern), ...
                                                    length(char(syncPattern)));
        end        
        
        function corrMatrix = calcCorrelationMatrix(obj, IQsamples, dimMatrix, ...
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

            delay_redFactor = fix(length(IQsamples) / dimMatrix);
            freq_redFactor = ceil((2 * maxDoppler / fresolution) / dimMatrix);
            Nsamples = delay_redFactor * dimMatrix;
            Nfrequencies = freq_redFactor * dimMatrix;
            
            obj.delays = (0:Nsamples-1)' / obj.fSampling;
            obj.dopplerFrequencies = linspace(-maxDoppler, maxDoppler, Nfrequencies);
            corrMatrix = zeros(dimMatrix);
            
            correlationMaximum = [0 0 0];
            %compute input fft
            input_fft = fft(IQsamples(:, 1) + 1i * IQsamples(:, 2), Nsamples);
            
            %cycle over frequencies
            for h = 1:Nfrequencies

                % DA RIVEDERE BISOGNA AGGIUNGERE ANCHE LA
                % MOLTIPLICAZIONE PER IL SYNC PATTERN, FORSE CONVIENE
                % FARE LA MOLTIPLICAZIONE PRIMA E POI ALLUNGARE O
                % RIDURRE IL TEMPO DI CHIP
                
                %TODO DT create correct interp PRN
                PRNsampled = obj.getPRNFromChipPeriods(shifts_nSamples_x_chipPeriod, segmentSize);
                PRN_Sync = PRNsampled.*obj.syncPattern;

                % taken from ambiguity function
                PRN_Sync_fft = fft(PRN_Sync.* ...
                                   exp(1i*2*pi*obj.dopplerFrequencies(h)* ...
                                       (obj.delays + currentSample / obj.fSampling)), ...
                                   Nsamples);
                product = PRN_Sync_fft .* conj(input_fft);
                columnCorrelation = abs(ifft(product));
                                
                % calculate and save peak precise position
                [max_columnCorrelation, maxIndex] = max(columnCorrelation, [], 1);
                if max_columnCorrelation > correlationMaximum(1)
                    correlationMaximum = [max_columnCorrelation maxIndex h];
                end
                
                %save reduced matrix, lineCorrelation is a column
                reducedLine = max(reshape(columnCorrelation, delay_redFactor, []),[], 1); %row vector
                corrMatrix(:, fix(h/freq_redFactor)) = max(corrMatrix(:, fix(h/freq_redFactor)), ...
                                                           reducedLine');
            end

            obj.fDoppler = obj.dopplerFrequencies(correlationMaximum(3));
            obj.startingSample = obj.delays(correlationMaximum(2)) * obj.fSampling + ...
                                 currentSample; % GIUSTO AGGIUNGERE IL CURRENT SAMPLE??
            obj.maxCorrelation = correlationMaximum(1);
        end

        function PRNsampled = getPRNFromChipPeriods(obj, nSamples_x_chipPeriod, windowSize)
            maxLength = nSamples_x_chipPeriod * length(obj.PRNsequence) * ...
                        obj.nPRN_x_Symbol * windowSize;
            PRNsampled = zeros(1, int32(maxLength + 0.5)); %ceil approx

            repSequence = repmat(obj.PRNsequence, 1, windowSize * obj.nPRN_x_Symbol);                
            PRNinterp = interp1(0:length(repSequence), [repSequence repSequence(1)], ...
                0:1 / nSamples_x_chipPeriod:length(repSequence), "previous"); %upsampling   
            repSize = size(PRNinterp,2);
            PRNsampled(period, 1:repSize) = PRNinterp;
        end
        
        %                           %
        %     DAVIDE TOMASELLA      %
        %  Functions for tracking.  %
        %  The parameter updates    %
        %  after each step.         %
        %                           %

        %DT$ already corrected
        function isAcquired = ifAcquiredFindCorrelationPeak(obj, corrMatrix, thresholdSTD, ...
                                                            outInterface)
            %ifAcquiredFindCorrelationPeak handles the decision if a peak
            %is found and the saving of the correlation peak parameters
            %   corrMatrix: correlation matrix (frequencies x time delays)
            %   thresholdSTD: n. STD over average to identify a correlation peak
            %   outInterface: outInterface.results contains the results

            % define dynamic threshold
            thresh = mean(corrMatrix, [1,2]) + thresholdSTD * std(corrMatrix, [1 2]);
            %find correlation maximum
            [max_corr, idx] = max(corrMatrix); %TODO use obj.maximumCorrelation
            if (max_corr > thresh)
                %get bet doppler and delay
                [idDopplerShift, idStartSample] = ind2sub(size(corrMatrix), idx);
                obj.fDoppler = obj.m_chipPeriods(idDopplerShift) - obj.txChipRate;
                obj.startingSample = obj.m_timeDelays(idStartSample) * obj.fSampling;
                %output saving
                outInterface.results.estimatedDoppler = obj.fDoppler;
                outInterface.results.estimatedDelay = obj.startingTime;
                %return successfull acquisition flag
                isAcquired = true;
            else
                %sync pattern not found not found
                isAcquired = false;
            end
        end
        
        %DT$ already corrected
        function obj = updateCorrelationPeak(obj, new_SamplesSymbolPeriod, advancement_startingSample)
            %updateCorrelationPeak handles the update of the peak position
            %during the tracking phase and keep updated the parameters
            %   new_SamplesSymbolPeriod: n. samples in new symbol period
            %   advancement_startingSample: n. samples to advance
            if new_SamplesSymbolPeriod<=0
                warning("Error. The number of samples per symbol period cannot be a negative value.")
            end
            if advancement_startingSample<=0
                warning("Error. The advancement size cannot be a negative value.")
            end
            %DT$ not needed thanks to dynamic properties
            %newChipPeriod = new_symbolPeriod / obj.nPRN_x_Symbol / obj.nChip_x_PRN;
            %obj.fDoppler = (1 / newChipPeriod) - obj.fModulation;
            obj.symbolPeriod = new_SamplesSymbolPeriod / obj.fSampling;
            obj.startingSample = obj.startingSample + uint32(advancement_startingSample);
        end

        %%%%%%%%% REAL PROPERTIES %%%%%%%%%

        %function set.fDoppler

        %function get.fDoppler 

        function set.startingSample(obj, iStartingSample)
            obj.startingSample = uint32(iStartingSample);
        end
       
        %function get.startingSample

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
