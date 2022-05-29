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
                figure(260);
                subplot(2,1,1);
                plot(reader.IQsamples_float(:,1));
                xlim([1 , length(reader.IQsamples)]);
                ylim([-1.1*max(reader.IQsamples(:,1)) 1.1*max(reader.IQsamples(:,1))]);
                title("Before down-conversion");
            end

            IQRef = obj.signalsCreation(obj.refAmplitude,obj.timebase(length(reader.IQsamples(:,1))),fdoppler,delay,phase);
            I = reader.IQsamples_float(:,1).*IQRef(:,1) - reader.IQsamples_float(:,2).*IQRef(:,2);
            Q = reader.IQsamples_float(:,1).*IQRef(:,2) + reader.IQsamples_float(:,2).*IQRef(:,1);
            reader.IQsamples_float = [I Q];
            
            if obj.DEBUG
                subplot(2,1,2);
                plot(reader.IQsamples_float(:,1));
                xlim([1 , length(reader.IQsamples)]);
                ylim([-1.1*max(reader.IQsamples(:,1)) 1.1*max(reader.IQsamples(:,1))]);
                title("After down-conversion");
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
                
                if obj.DEBUG
                    figure(250);
                    plot(reader.IQsamples_float(:,1));
                    hold on;
                end

                fNyq = obj.fsampling / 2;
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
%                 [ord , W] = cheb1ord(passBand/fNyq,stopBand/fNyq,obj.ripple,stopBand_attenuation);

                myDEBUG = false;
                if myDEBUG
                    disp("Attenuazione " + num2str(stopBand_attenuation)); %% PER DEBUG
                    disp("Ordine " + num2str(ord) + " Banda " + num2str(W)); %% PER DEBUG
                    b = fir1(ord,W,"low");
                    freqz(b);
                end
                
                [b , a] = butter(ord,W);
%                 [b , a] = cheby1(ord,obj.ripple,W);
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
                    plot(reader.IQsamples_float(:,1));
                    hold off;
                    xlim([length(reader.IQsamples)/2 , length(reader.IQsamples)/2+500]);
                    title("Signal filtering");
                    legend("Before","After");
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