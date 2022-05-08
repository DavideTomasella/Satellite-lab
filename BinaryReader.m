%
% First implementation: Mattia Piana
% Review and Testing: Davide Tomasella
%
classdef BinaryReader < handle
    %BinaryReader handles...
   properties (Access = { ?classUnderTest, ?matlab.unittest.TestCase })
        isReadFileConfigured=false;
    end
   properties (SetAccess=private, GetAccess=public)
       fullfilename
       nByte_per_sample uint32{mustBeGreaterThan(nByte_per_sample,0),mustBeLessThanOrEqual(nByte_per_sample,4)}
       nSamples
       IQsamples
       IQsamples_float
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
            if nargin < 4
                nBit_per_sample = 16;
            end
            if nargin < 5
                relativePath = true;
            end
            if relativePath
                dirName = fullfile(".",dirName);
            end
            obj.fullfilename = fullfile(dirName,fileName);
            %DT$ approx at closest multiple of 8 -> Matlab limitation!
            obj.nByte_per_sample = int16(nBit_per_sample/8);
            obj.nSamples = obj.getNSamplesFromFile();
            obj.isReadFileConfigured = true;
        end

        %                        %
        % ACTION METHODS         %
        %                        %

        function n = getNSamplesFromFile(obj)
            %DT$ ftell return the dimension in bytes + "/2"= 2 components
            %fseek(fid, 0, 'eof');
            %n = ftell(fid)/obj.nByte_per_sample/2;
            try
                s = dir(obj.fullfilename);         
                n = s.bytes/obj.nByte_per_sample/2;
            catch
                disp("Error file not found: aborted");
             end
        end

        function obj = readFile(obj,currentSample,nSamples)
            %Function readFile that reads a window of the binary file spanning from currentSample
            %to currentSample+nSamples and instantiate the field IQsamples
            %with the int16 data
            if currentSample+nSamples>obj.nSamples
                currentSample=obj.nSamples-nSamples;
            end
            if currentSample < 0
                disp("Error in file reading starting point: set to 0")
                currentSample = 0;
            end
            tmpData = obj.openReadCloseBinaryFile(currentSample,nSamples);
            %obj.IQsamples = obj.formatSamples(tmpData);
            obj.IQsamples=tmpData;
            clear tmpData
        end

        function oSamples_float = get.IQsamples_float(obj)
            oSamples_float = single(obj.IQsamples);
        end

        %                        %
        % PRIVATE METHODS        %
        %                        %
     
        function tmpData = openReadCloseBinaryFile(obj,currentSample,nSamples)
            try
                fid = fopen(obj.fullfilename);
                if(fid~=-1)
                    fseek(fid, currentSample*obj.nByte_per_sample*2, 'bof');
                    tmpData = fread(fid,[2 nSamples],obj.getSampleFormat());
                else 
                    disp("Matlab cannot open the file: aborted",obj.fullfilename);
                end
                fclose(fid);
            catch
                disp("Error in opening or closing the file: aborted");
            end
        end

        function format = getSampleFormat(obj)
            if obj.nByte_per_sample == 1
                format = "int8=>int16";
            elseif obj.nByte_per_sample == 2
                format = "int16=>int16";
            elseif obj.nByte_per_sample == 3
                format = "bit24=>int32";
            elseif obj.nByte_per_sample == 4
                format = "int32=>int32";
            else
                disp("Errore: impossible value for nByte_per_sample")
            end
        end
   end
end