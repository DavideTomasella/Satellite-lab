%
% First implementation: Gabriele Andreetta
% Review and Testing: Gabriele Andreetta
%
classdef DownconverterFilter < handle
    %DownconverterFilter handles...
    

    properties (SetAccess=private, GetAccess=public)
        fsampling
        ripple
        attenuation
        refAmplitude
    end

    properties (SetAccess=public, GetAccess=public)
        DEBUG
    end

    methods
        function obj = DownconverterFilter(DEBUG)
            %DownconverterFilter constructor
            %   Create DownconverterFilter class
            if nargin < 1
                DEBUG = false;
            end
            obj.DEBUG = DEBUG;
        end

        function obj = configDownConverter(obj,fSampling)
            obj.fsampling = fSampling;
            obj.refAmplitude = 1;
        end
        
        % Down convertion achieved by multiplication with the complex exponential
        % exp(-1i*2*pi*fdoppler*(t+delay))
        function reader = downConverter(obj,reader,fdoppler,delay,phase)
            if obj.DEBUG
                dcf1 = figure(250);
                movegui(dcf1,"northwest")
                plot(obj.timebase(length(reader.IQsamples(:,1))),reader.IQsamples_float(:,1));
                grid on;
                xlim tight;
                ylim([-1.1*max(reader.IQsamples(:,1)) 1.1*max(reader.IQsamples(:,1))]);
                xlabel("Time [s]");
                ylabel("Amplitude");
                title("Before down-conversion");
%                 savePdf(dcf1,"DownConverter1",true);
            end

            IQRef = obj.signalsCreation(obj.refAmplitude,obj.timebase(length(reader.IQsamples(:,1))),fdoppler,delay,phase);
            I = reader.IQsamples_float(:,1).*IQRef(:,1) - reader.IQsamples_float(:,2).*IQRef(:,2);
            Q = reader.IQsamples_float(:,1).*IQRef(:,2) + reader.IQsamples_float(:,2).*IQRef(:,1);
            reader.IQsamples_float = [I Q];
            
            if obj.DEBUG
                dcf2 = figure(251);
                movegui(dcf2,"northwest")
                plot(obj.timebase(length(reader.IQsamples(:,1))),reader.IQsamples_float(:,1));
                grid on;
                xlim tight;
                ylim([-1.1*max(reader.IQsamples(:,1)) 1.1*max(reader.IQsamples(:,1))]);
                xlabel("Time [s]");
                ylabel("Amplitude");
                title("After down-conversion");
%                 savePdf(dcf1,"DownConverter2",true);
                pause(0.3)
            end

            clear IQRef
            clear I
            clear Q
        end

        function obj = configFilter(obj,ripple_dB,attenuation_dB_dec,fSampling)
            if nargin < 4 && isempty(obj.fsampling)
                disp("Error: sample frequency not present");
            elseif nargin == 4
                obj.fsampling = fSampling;
            end
            if ripple_dB <= 0
                disp("Error: ripple must be higher than 0, default value 1dB");
                obj.ripple=1;
            else
                obj.ripple = ripple_dB;
            end   
            if attenuation_dB_dec <= 0
                disp("Error: stopband attenuation must be higher than 0, default value 20dB");
                obj.attenuation = 20;
            else
                obj.attenuation = attenuation_dB_dec;
            end
        end

        function reader = downFilter(obj,reader,passBand,chipFrequency)
            if passBand > 0
                
                fNyq = obj.fsampling / 2;
                if obj.DEBUG
                    ff1 = figure(252);
                    movegui(ff1,"southwest")
                    plot(obj.timebase(length(reader.IQsamples(:,1))),reader.IQsamples_float(:,1));
                    grid on;
                    xlim([length(reader.IQsamples)/2/obj.fsampling , length(reader.IQsamples)/2/obj.fsampling+500/obj.fsampling]);
                    if true
                        ff2 = figure(253);
                        movegui(ff2,"south")
                        plot((0:length(reader.IQsamples(:,1))-1)*obj.fsampling/length(reader.IQsamples(:,1))-fNyq,abs(fftshift(fft(reader.IQsamples_float(:,1)))));
                        grid on;
                        xlim tight;
                    end
                end

                if passBand >= fNyq - 1
                    passBand = (fNyq - 1) / 2;
                    disp("Error: passband overcomes the Nyquist frequency fs/2, used default passband: "+num2str(passBand,5));
                end
                stopBand = fNyq-1;
                stopBand_attenuation = obj.attenuation * log10(stopBand / passBand);
                
                if stopBand_attenuation > 3000 %dB
                    if obj.DEBUG
                        disp("Warning: stopband attenuation is too high, dark band values already null. " + ...
                             "Impossible calculate filter order. stopBand at 3000dB");
                    end
                    stopBand_attenuation = 3000;
                end

                [ord , W] = buttord(passBand/fNyq,stopBand/fNyq,obj.ripple,stopBand_attenuation);
                % [ord , W] = cheb1ord(passBand/fNyq,stopBand/fNyq,obj.ripple,stopBand_attenuation);
                
                [b , a] = butter(ord,W);
                % [b , a] = cheby1(ord,obj.ripple,W);
                [d , w] = grpdelay(b,a);
                IQfiltered = filter(b,a,reader.IQsamples_float,[],1);
    
                if chipFrequency >= fNyq
                    disp("Warning: undersampling");
                end

                delay_index = find(w <= chipFrequency/fNyq , 1 , "last");
                samples_delay = round(d(delay_index));
                IQfiltered = [IQfiltered(samples_delay+1:end,:) ; zeros(samples_delay,2)];
                reader.IQsamples_float = IQfiltered;

                if obj.DEBUG
                    figure(ff1);
                    hold on;
                    plot(obj.timebase(length(reader.IQsamples(:,1))),reader.IQsamples_float(:,1));
                    hold off;
                    xlabel("Time [s]");
                    ylabel("Amplitude");
                    title("Filter's effect on time domain signal");
                    legend("Before filtering","After filtering");
                    if true
                        figure(ff2);
                        hold on;
                        plot((0:length(reader.IQsamples(:,1))-1)*obj.fsampling/length(reader.IQsamples(:,1))-fNyq,abs(fftshift(fft(reader.IQsamples_float(:,1)))));
                        hold off;
                        xlabel("Frequency [Hz]");
                        ylabel("Amplitude");
                        title("Filter's effect on signal spectrum");
                        legend("Before filtering","After filtering");
%                         savePdf(ff1,"Filter");
%                         savePdf(ff2,"FilterFFT");
                    end
                    pause(0.3)
                end

                clear IQfiltered
                clear stopBand
                clear stopBand_attenuation
                clear ord
                clear W
                clear fNyq
                clear delay_index
                clear samples_delay
                clear b
                clear d
                clear w
            else
                if obj.DEBUG
                    disp("Filter absent");
                end
            end
        end
    end

    methods(Access = private)
        function time = timebase(obj,Nsamples)
            time = (0:Nsamples-1)'/obj.fsampling;
        end

        function IQRef = signalsCreation(~,amplitude,time,fdoppler,delay,phase)
            ampl = double(amplitude);
            IRef = ampl*cos(2*pi*fdoppler*(time+delay)+phase);
            QRef = -ampl*sin(2*pi*fdoppler*(time+delay)+phase);
            IQRef = [IRef , QRef];
            clear ampl
            clear IRef
            clear QRef
        end
    end
end