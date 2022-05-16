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
            Nsamples = size(IQsamples,1);
            obj.m_chipPeriods = obj.txChipRate - maxDoppler:fresolution:obj.txChipRate + maxDoppler;
            obj.m_timeDelays = (currentSample + (1:Nsamples))' / obj.fSampling; %column vector
            
            %DT$ error see upsampling_examples.m
            PRNsampled = repmat(obj.PRNsequence, 1, obj.nPRN_x_Symbol);
            
            % note: il chip period varia in funzione del doppler che stai
            % considerando attraverso il nuovo chiprate (i.e. modulation
            % frequency of PRN bits) saved in obj.frequencies

            % note: la corelazione va fatta anche sul pattern di sync!!!
            
            % nel caso invece volessi utilizzare il symbol period ti ho
            % creato la prperty obj.txSymbolRate = chipRate / nChip_x_PRN / nPRN_x_Symbol;
            % Modifica quindi obj.frequencies
            % Ricordati però di modificare l'output in ifAcquiredFindCorrelationPeak
            
            % Per decidere cosa ti serve devi capire cosa è questo exp(1i*2*pi*f*time)
            % magari in realtà qui va solo il doppler shift e non la total
            % frequenza di modulation del nostro segnale

            corrMatrix = zeros(Nsamples, length(obj.m_chipPeriods));
            j=1;
            for f = obj.m_chipPeriods
                x_t = fft((IQsamples(1,:)+1i*IQsamples(2,:)).*exp(1i*2*pi*f*time));
                x_tau = conj(fft(PRNsampled));
                corrMatrix(:,j)=fftshift(abs(ifft(x_t.*x_tau)));
                j=j+1;
            end
            %note: test the algorithm with real signals
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
