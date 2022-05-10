%
% First implementation: Marcello Raimondi
% Review and Testing:
%
classdef Correlator < handle
    %Correlator handles...
    
    properties
        
    end
        
    methods
        function obj = Correlator()
            f_linspace
            PRN
        end
        
        function obj = configCorrelatorMatrix()
            % setup the indeces for the matrix and
            % import the PRN sequences from inout
        end
        
        function corrMatrix = calcCorrelationMatrix(samples)
            
            N = size(samples[1]);
            matr = zeros(1, 2*N-1)
            time = linspace(1, N, N);
            j=1;
            for f = obj.f_linspace
                x_t = fft(samples[0].*exp(i*2*pi*f*time)+samples[1].*exp((i*2*pi*f*time)+pi/2));
                x_tau = conj(fft(obj.PRN));
                    
                corr_vec = xcorr(x_t, x_tau);
                matr(j)=corr_vec;
            end
            
        end
        
        function obj = findCorrelationPeak(matr, samples, inout)
            [~, idx] = max(matr);
        end
end
