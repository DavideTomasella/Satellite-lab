%
% First implementation: Mattia Piana
% Review and Testing:
%
classdef BinaryReader < handle
    %BinaryReader handles...
    properties (Access = { ?classUnderTest, ?matlab.unittest.TestCase })
        isReadFileConfigured=false;
    end
   properties (SetAccess=private, GetAccess=public)
       fullfilename
       nBit_per_sample
       nSamples
       nSamples_float
       IQsamples
   end
    
   methods
       function obj = BinaryReader()
       end

        %                        %
        % CONFIGURATION METHODS  %
        %                        %

        function obj = configReadFile(obj,dirName,fileName,nBit_per_sample,relativePath)
            %configCreateSettings transform the json input file into a config
            %struct and validate the input parameters
            %   dirName: string with absolute or relative path where
            %   configuration file and sampled data are stored
            if nargin < 5
                relativePath = true;
            end
            if relativePath
                dirName = fullfile(".",dirName);
            end
            obj.fullfilename = fullfile(dirName,fileName);
            obj.nBit_per_sample = nBit_per_sample;
            obj.nSamples = obj.getNSamplesFromFile(obj.fullfilename);
            obj.isReadFileConfigured = true;
        end

        function obj = readFile(obj,currentSample,nSamples)
            if currentSample+nSamples>obj.nSamples
                currentSample=obj.nSamples-nSamples;
            end
            if currentSample < 0
                disp("erroooor")
                return
            end
            tmpData = obj.openReadCloseBinaryFile(currentSample*2*obj.nBit_per_sample,nSamples*2*obj.nBit_per_sample);
            obj.IQsamples = obj.formatSamples(tmpData);
            clear tmpData
        end
        %implement getNSamplesFromFile (dimensione... numero bits non è
        %divisibile... )
        %check currentSample*2*obj.nBit_per_sample
        %implement openReadCloseBinaryFile (try catch... verificare file
        %corretto... ) -> fopen, fseek+read(length)/read(skip,length), close
        %implement (maybe) formatSamples (remove average...)
        %method get.nSamples_float: ritorna nSamples convertito in single
   end
end