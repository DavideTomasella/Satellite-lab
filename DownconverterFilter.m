%
% First implementation: Gabriele Andreetta
% Review and Testing:
%
classdef DownconverterFilter < handle
    %DownconverterFilter handles...
    

    properties (SetAccess=public, GetAccess=public)
        fsampling
        ripple
        passbandstopbandratio
        attenuation
    end

    methods
        function obj = DownconverterFilter()
            %DownconverterFilter constructor
            %   Create DownconverterFilter class
        end

        function obj = configDownConverter(obj,fSampling)
            obj.fsampling = fSampling;
        end
        
        function time = timebase(obj,Nsamples)
            time = (0:Nsamples-1)'/obj.fsampling;
        end

        function IQRef = signalsCreation(~,time,fdoppler,delay)
            IRef = cos(2*pi*fdoppler*(time+delay));
            QRef = -1*sin(2*pi*fdoppler*(time+delay));
            IQRef = [IRef , QRef];
            clear IRef
            clear QRef
        end
        
        % Down convertion achieved by multiplication with the complex exponential
        % exp(-1i*2*pi*fdoppler*(t+delay))
        function downSamples = downConverter(obj,IQsamples,fdoppler,delay)
            IQRef = obj.signalsCreation(obj.timebase(length(IQsamples(:,1))),fdoppler,delay);
            I = IQsamples(:,1).*IQRef(:,1) - IQsamples(:,2).*IQRef(:,2);
            Q = IQsamples(:,1).*IQRef(:,2) + IQsamples(:,2).*IQRef(:,1);
            downSamples = [I Q];
            clear IQRef
            clear I
            clear Q
        end

        function obj = configFilter(obj,passbandstopbandratio,ripple_dB,attenuation_dB)
            if ripple_dB <= 0
                disp("Error: ripple must be different from 0, default value 1dB");
                obj.ripple=1;
            else
                obj.ripple = ripple_dB;
            end            
            obj.passbandstopbandratio = passbandstopbandratio;
            obj.attenuation = attenuation_dB;
        end

        function filteredSamples = downFilter(obj,IQ,passBand)
            stopband = obj.passbandstopbandratio * passBand;
            if stopband > obj.fsampling / 2
                stopband = obj.fsampling / 2;
                passBand = stopband / obj.passbandstopbandratio;
                disp("Error: stopband overcomes the Nyquist frequency fs/2, used max bandwidth: "+num2str(passBand,5));
            end
            h = fdesign.lowpass(passBand, stopband, obj.ripple, obj.attenuation, obj.fsampling);
            Hd = design(h, 'butter', 'MatchExactly', 'passband');   % match the passband frequency
            Ifiltered = filter(Hd,IQ(:,1));
            Qfiltered = filter(Hd,IQ(:,2));
            filteredSamples = [Ifiltered, Qfiltered];
            clear Ifiltered
            clear Qfiltered
        end
    end
end