%
% First implementation: Gabriele Andreetta
% Review and Testing:
%
classdef DownconverterFilter < handle
    %DownconverterFilter handles...
    

    properties (SetAccess=private, GetAccess=public)
        fsampling
        ripple
        passbandstopbandratio
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

        function obj = configFilter(obj,passbandstopbandratio,ripple_dB,attenuation_dB,fSampling)
            if nargin < 6 && isempty(obj.fsampling)
                disp("Error: sample frequency not present");
            end
            if ripple_dB <= 0
                disp("Error: ripple must be higher than 0, default value 1dB");
                obj.ripple=1;
            else
                obj.ripple = ripple_dB;
            end   
            if passbandstopbandratio <= 1
                disp("Error: passband-stopband ratio must be higher than 1, default value 10");
                obj.passbandstopbandratio = 10;
            else
                obj.passbandstopbandratio = passbandstopbandratio;
            end
            if attenuation_dB <= 0
                disp("Error: stopband attenuation must be higher than 0, default value 20dB");
                obj.attenuation = 20;
            else
                obj.attenuation = attenuation_dB;
            end
            if attenuation_dB > 3000
                disp("Error: stopband attenuation is too high, max value 3000dB");
                obj.attenuation = 3000;
            else
                obj.attenuation = attenuation_dB;
            end
        end

        function reader = downFilter(obj,reader,passBand,chipFrequency)
            fNyq = obj.fsampling / 2;
            if passBand >= fNyq - 1
                passBand = (fNyq - 1) / obj.passbandstopbandratio;
                disp("Error: passband overcomes the Nyquist frequency fs/2, used default passband: "+num2str(passband,5));
            end
            stopband = obj.passbandstopbandratio * passBand;
            if stopband > fNyq
                stopband = fNyq-1;
%                 passBand = stopband / obj.passbandstopbandratio;
                disp("Error: stopband overcomes the Nyquist frequency fs/2, used max stopband: "+num2str(stopband,5));
            end
            [ord , W] = buttord(passBand/fNyq,stopband/fNyq,obj.ripple,obj.attenuation);
            [b , a] = butter(ord,W);
            [d , w] = grpdelay(b,a);

            Ifiltered = filter(b,a,reader.IQsamples_float(:,1)');
            Qfiltered = filter(b,a,reader.IQsamples_float(:,2)');
            
            if chipFrequency >= fNyq
                disp("Warning: undersampling");
            end
            delay_index = find(w <= chipFrequency/fNyq);
%             fchip_index = delay_index(end);
            samples_delay = int16(d(delay_index(end)));
            Ifiltered = circshift(Ifiltered,-samples_delay);
            Qfiltered = circshift(Qfiltered,-samples_delay);

%             h = fdesign.lowpass(passBand, stopband, obj.ripple, obj.attenuation, obj.fsampling);
%             Hd = design(h, 'butter', 'MatchExactly', 'passband');   % match the passband frequency
%             Ifiltered = filter(Hd,reader.IQsamples(:,1));
%             Qfiltered = filter(Hd,reader.IQsamples(:,2));

            reader.IQsamples = int16([Ifiltered', Qfiltered']);
            clear Ifiltered
            clear Qfiltered
            clear fNyq
        end
    end
end