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
            IRef = sin(2*pi*fdoppler*(time+delay));
            QRef = cos(2*pi*fdoppler*(time+delay));
            IQRef = [IRef , QRef];
            clear IRef
            clear QRef
        end
        
        % Down convertion achieved by multiplication with the carrier at
        % the doppler frequency and the time delay.
        % This cannot be used with a filter because of the aliasing between
        % the positive and negative spectrum, a filter would remove the
        % code, not the carrier
        function downSamples = downConverter(obj,IQsamples,fdoppler,delay)
            downSamples = IQsamples.*obj.signalsCreation(obj.timebase(length(IQsamples(:,1))),fdoppler,delay);
        end

        function obj = configFilter(obj,passbandstopbandratio,ripple_dB,attenuation_dB)
            obj.ripple = ripple_dB;
            obj.passbandstopbandratio = passbandstopbandratio;
            obj.attenuation = attenuation_dB;
        end

        function filteredSamples = downFilter(obj,IQ,passBand)
            stopband = obj.passbandstopbandratio*passBand;
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