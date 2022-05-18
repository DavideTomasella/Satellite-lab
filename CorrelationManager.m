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
        %matrix dimensions
        m_chipPeriods
        m_timeDelays

        dopplerFrequencies
        delays;
    end
    properties (SetAccess=public, GetAccess=public)
        %output real properties
        fDoppler {single}
        startingSample {uint32}
        maxCorrelation
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
        
        function corrMatrix = calcCorrelationMatrix(obj, IQsamples, ...
                                                    maxDoppler, fresolution, currentSample)
            %calcCorrelationMatrix handles the creation of the correlation
            %matrix (doppler frequencies and time delays)
            % ##Code adapted from:
            % ##source: http://maxwell.sze.hu/~ungert/Radiorendszerek_satlab/Segedanyagok/Ajanlott_irodalom/matlab_simulations_for_radar_systems_design_(cuppy).pdf
            % ##page: 247/686
            %   IQsamples: column vector with input samples
            %   maxDoppler: maximum doppler
            %   fresolution: doppler frequency step size
            %   currentSample: delay (in samples) of the considered signal
            
            obj.dopplerFrequencies = -maxDoppler:fresolution:maxDoppler;
            Nfrequencies = length(obj.dopplerFrequencies);
            Nsamples = length(IQsamples);
            obj.delays = (0:Nsamples-1)/obj.fSampling;
            
            frequency_reductionFactor = 10;
            delay_reductionFactor = 10;
            corrMatrix = zeros(Nsamples/delay_reductionFactor,length(obj.dopplerFrequencies)/frequency_reductionFactor);
            partialMatrix = zeros(Nsamples/delay_reductionFactor,frequency_reductionFactor);
            
            max_single_correlation = 0;
            correlationMaximum = [0 0 0];
            
            for h = 1:(Nfrequencies/frequency_reductionFactor)
                for j=1:frequency_reductionFactor
                    frequency_index = j+frequency_reductionFactor*(h-1);
                    
                    % DA RIVEDERE BISOGNA AGGIUNGERE ANCHE LA
                    % MOLTIPLICAZIONE PER IL SYNC PATTERN, FORSE CONVIENE
                    % FARE LA MOLTIPLICAZIONE PRIMA E POI ALLUNGARE O
                    % RIDURRE IL TEMPO DI CHIP
                    
                    PRNsampled = obj.getPRNFromChipPeriods(shifts_nSamples_x_chipPeriod, segmentSize);
                    PRN_Sync = PRNsampled.*obj.syncPattern;
    %                 x_t = fft((IQsamples(1,:)+1i*IQsamples(2,:)).*exp(1i*2*pi*f*time));
    %                 x_tau = conj(fft(PRNsampled));

                    % taken from ambiguity function
                    PRN_Sync_fft = fft(PRN_Sync);   %aggiungere il PADDING
                    input_fft = fft(IQsamples.*exp(1i*2*pi*obj.dopplerFrequencies(frequency_index)*delay));
                    product = input_fft .* conj(PRN_Sync_fft);
    
                    single_correlation=fftshift(abs(ifft(product)));
                    [max_single_correlation, maxindex] = max(single_correlation);

                    if max_single_correlation > correlationMaximum(1)
                        correlationMaximum = [max_single_correlation maxindex frequency_index];
                    end
                    
                    for k = 1:(Nsamples/delay_reductionFactor)
                        index1 = (k-1)*delay_reductionFactor+1;
                        index2 = k*delay_reductionFactor;
                        partialMatrix(k,j) = max(single_correlation(index1:index2));
                    end

%                     correlazione;reshape(correlazione,10,[]);max(correlazione,1);
%                     maxdownsampled=time/10;
%                     corrMatrix(:,fix(j/10))=max(maxdownsampled,corrMatrix(:,fix(j/10)));
                end
                corrMatrix(:,h) = max(partialMatrix,[],2);  % column vector with maximum for each row
            end
            obj.fDoppler = obj.dopplerFrequencies(correlationMaximum(3));
            obj.startingSample = correlationMaximum(2) + currentSample; % GIUSTO AGGIUNGERE IL CURRENT SAMPLE??
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
            [max_corr, idx] = max(corrMatrix);
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
