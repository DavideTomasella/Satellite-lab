%
% First implementation: Gabriele Andreetta
% Review and Testing: Gabriele Andreetta
%
classdef DownconverterFilter < handle
    %DownconverterFilter handles...
    

    properties (SetAccess=private, GetAccess=private)
        fsampling
        ripple
        attenuation
        refAmplitude
    end

    methods
        function obj = DownconverterFilter()
            %DownconverterFilter constructor
            %   Create DownconverterFilter class
        end

        function obj = configDownConverter(obj,fSampling)
            obj.fsampling = fSampling;
            obj.refAmplitude = 1;
        end
        
        % Down convertion achieved by multiplication with the complex exponential
        % exp(-1i*2*pi*fdoppler*(t+delay))
        function reader = downConverter(obj,reader,fdoppler,delay)
            IQRef = obj.signalsCreation(obj.refAmplitude,obj.timebase(length(reader.IQsamples(:,1))),fdoppler,delay);
            I = reader.IQsamples_float(:,1).*IQRef(:,1) - reader.IQsamples_float(:,2).*IQRef(:,2);
            Q = reader.IQsamples_float(:,1).*IQRef(:,2) + reader.IQsamples_float(:,2).*IQRef(:,1);
            reader.IQsamples_float = [I Q];
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
            fNyq = obj.fsampling / 2;
            if passBand >= fNyq - 1
                passBand = (fNyq - 1) / 2;
                disp("Error: passband overcomes the Nyquist frequency fs/2, used default passband: "+num2str(passBand,5));
            end
            stopBand = fNyq-1;
            stopBand_attenuation = obj.attenuation * stopBand/passBand;
            if stopBand_attenuation > 3000
                disp("Error: stopband attenuation is too high, max value 3000dB");
                stopBand_attenuation = 3000;
            end
            [ord , W] = buttord(passBand/fNyq,stopBand/fNyq,obj.ripple,stopBand_attenuation);
            b = fir1(ord,W,"low");
%             freqz(b);
            [d , w] = grpdelay(b);
            IQfiltered = filter(b,1,reader.IQsamples_float,[],1);

            if chipFrequency >= fNyq
                disp("Warning: undersampling");
            end
            delay_index = find(w <= chipFrequency/fNyq , 1 , "last");
            samples_delay = round(d(delay_index));
            IQfiltered = [IQfiltered(samples_delay+1:end,:) ; zeros(samples_delay,2)];
            reader.IQsamples_float = IQfiltered;
            clear IQfiltered
            clear stopBand
            clear ord
            clear W
            clear fNyq
            clear delay_index
            clear samples_delay
            clear b
            clear d
            clear w
        end
    end

    methods(Access = private)
        function time = timebase(obj,Nsamples)
            time = (0:Nsamples-1)'/obj.fsampling;
        end

        function IQRef = signalsCreation(~,amplitude,time,fdoppler,delay)
            ampl = double(amplitude);
            IRef = ampl*cos(2*pi*fdoppler*(time+delay));
            QRef = -ampl*sin(2*pi*fdoppler*(time+delay));
            IQRef = [IRef , QRef];
            clear ampl
            clear IRef
            clear QRef
        end
    end
end