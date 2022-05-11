%
% First implementation: Marcello Raimondi
% Review and Testing:
%
classdef Correlator < handle
    %Correlator handles...
    
    properties
        fSampling
        nPRN_x_Symbol
        PRNcode
        fModulation
        frequencies
        delays        
    end
        
    methods
        function obj = Correlator()
            
        end        
        
        function obj = configCorrelatorMatrix(obj,fSampling,nPRN_x_Symbol,PRNcode,fModulation)
            obj.fSampling = fSampling;
            obj.nPRN_x_Symbol = nPRN_x_Symbol;
            obj.PRNcode = PRNcode;
            obj.fModulation =fModulation;
        end        
        
        function corrMatrix = calcCorrelationMatrix(obj,IQsamples,maxDoppler,fresolution)
            
            Nsamples = size(IQsamples,1);
            obj.frequencies = obj.fModulation-maxDoppler:fresolution:obj.fModulation+maxDoppler;
            corrMatrix = zeros(Nsamples, length(obj.frequencies));
            obj.delays = (1:Nsamples)'/obj.fSampling;
            %DT$ TODO create PRNsequence (time-continuous vector) from obj.PRNcode (bit sequence)!
            j=1;
            for f = obj.frequencies
                x_t = fft((IQsamples(1,:)+1i*IQsamples(2,:)).*exp(1i*2*pi*f*time));
                x_tau = conj(fft(PRNsequence));
                corrMatrix(:,j)=fftshift(abs(ifft(x_t.*x_tau)));
                j=j+1;
            end
        end
        
        function isAcquired = ifAcquiredFindCorrelationPeak(obj,corrMatrix,thresholdSTD,currentSample,oResults)

            thresh = mean(corrMatrix, 'all') + thresholdSTD * std(corrMatrix, 'all');
            [max_corr, idx] = max(corrMatrix);
            if (max_corr > thresh)
                [idDopplerShift,idStartSample] = ind2sub(size(corrMatrix),idx);
                oResults.estimatedDoppler = obj.frequencies(idDopplerShift);
                oResults.estimatedDelay = oResults.estimatedDelay + ...
                                          (currentSample + obj.delays(idStartSample))*obj.fSampling;
                isAcquired = true;
            else
                isAcquired = false;
            end
        end
    end
end
