%
% First implementation: Marcello Raimondi
% Review and Testing:
%
classdef Correlator < handle
    %Correlator handles...
    
    properties
        sampling_freq
        PRN
        t_delay
        dopp_freq
        f_linspace
    end
        
    methods
        function obj = Correlator()
            
        end
        
        
        function obj = configCorrelatorMatrix()
            sampling_freq = 10e6
            PRN = [0 1 0]   %placeholder
            f_linspace = 0.823e6:50e3:1.223e6
        end
        
        
        function matr = calcCorrelationMatrix(samples)
            
            N = size(samples[1]);
            matr = zeros(1, 2*N-1);
            time = linspace(1, N, N);
            j=1;
            for f = obj.f_linspace
                x_t = fft(samples[1,:].*exp(i*2*pi*f*time)+samples[2,:].*exp((i*2*pi*f*time)+pi/2));
                x_tau = conj(fft(obj.PRN));
                    
                matr(j)=ifft(xcorr(x_t, x_tau));
                j++;
            end
            
        end
        
        
        function obj = findCorrelationPeak(matr, samples, inout)
            inout.estimatedDoppler = dopp_freq;
            inout.estimatedDelay = t_delay;
        end
        
        
        function acq = isSignalAcquired(matr)
            [max_corr, idx] = max(matr);
            if(max_corr >= obj.thresh)
                dopp_freq = idx[1];
                t_delay = idx[2]*sampling_freq;
                acq = true
            else
                acq = false
        end
end
