%
% First implementation: Gabriele Andreetta
% Review and Testing:
%
classdef DownconverterFilter < handle
    %DownconverterFilter handles...
    

    properties (SetAccess=public, GetAccess=public)
        fsampling
        % Needed for filter
        lowPassFilter_f
        isLowPassFilterConfigured
    end

    methods
        function obj = DownconverterFilter()
            %DownconverterFilter constructor
            %   Create DownconverterFilter class
        end

        function obj = configDownConverter(obj,fsampling)
            obj.fsampling = fsampling;
        end
        
        function time = timebase(obj,Nsamples)
            time = (0:Nsamples)/obj.fsampling;
            time = time';
        end

        function IQRef = signalsCreation(obj,time,fdoppler,delay)
            IRef = sin(2*pi*fdoppler*(time+delay));
            QRef = cos(2*pi*fdoppler*(time+delay));
            IQRef = [IRef , QRef];
        end

        function obj = downConverter(obj,IQsamples,fdoppler,delay)
            obj.downconvertedI = IQsamples.*signalsCreation(timebase(length(IQsamples(:,1))),fdoppler,delay);
        end

        function filteredSamples = filter(downIQ,passBand)
            
        end
        
    end
end