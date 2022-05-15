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
       binDirectory
       inFullfilename
       nByte_per_sample uint32{mustBeGreaterThan(nByte_per_sample,0),mustBeLessThanOrEqual(nByte_per_sample,4)}
       totSamples
   end
   properties (SetAccess=public, GetAccess=public)
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
            obj.binDirectory = dir(dirName);
            obj.inFullfilename = fullfile(dirName,fileName);
            %DT$ approx at closest multiple of 8 -> Matlab limitation!
            obj.nByte_per_sample = int16(nBit_per_sample/8);
            obj.totSamples = obj.getNSamplesFromFile();
            obj.isReadFileConfigured = true;
        end

        %                        %
        % ACTION METHODS         %
        %                        %

        function n = getNSamplesFromFile(obj)
            fid = fopen(obj.inFullfilename);
            try
                if(fid ~= -1)
                    fseek(fid, 0, 'eof');
                    n = ftell(fid) / 2 / obj.nByte_per_sample;
                    if fclose(fid) == -1
                        disp("Error: impossible close file")
                    end
                else
                   disp("Error in opening the file");
                end
            catch
                warning("Error: impossible read file %s",obj.inFullfilename);
            end
        end

        function obj = readFile(obj, currentSample, nSamples)
            %Function readFile that reads a window of the binary file spanning from currentSample
            %to currentSample+nSamples and instantiate the field IQsamples
            %with the int16 data
            if currentSample + nSamples > obj.totSamples
                disp("Not enough samples, padding required");
            end
            if currentSample < 0
                disp("Error in file reading starting point: set to 0")
                currentSample = 0;
            end
            toReadData = min(nSamples, obj.totSamples - currentSample);
            tmpData = obj.openReadCloseBinaryFile(currentSample, toReadData);
            padding = zeros(2, nSamples - toReadData);
            %obj.IQsamples = obj.formatSamples(tmpData);
            obj.IQsamples = [tmpData padding]';
            clear tmpData
        end

        function fullname = saveToBynaryFile(obj,samples,filename,append)
        %Append at the end of a Binary File specified by filename the
        %content of the n-bit data in the coloumn vector samples returns an
        %log error if the saving failed
            try
                folder = obj.binDirectory.folder;
                fullname = fullfile(folder,filename);
                if nargin<4 
                    append = false;
                end
                if append
                    fileID = fopen(fullname,'a');
                else
                    fileID = fopen(fullname,'w');
                end
                fwrite(fileID,samples',obj.getSavingFormat());
                if fclose(fileID)==-1
                    disp("Error: impossible close file")
                end
            catch
                warning("Error: impossible save file %s",fullname);
            end
        end

        function oSamples_float = get.IQsamples_float(obj)
            oSamples_float = single(obj.IQsamples);
        end

        function set.IQsamples_float(obj, iSamples_float)
            %corretto, IQsamples_float is a virual property
            % the cast grants the convertion tu the correct types
            obj.IQsamples = cast(iSamples_float,obj.getSavingFormat());
        end

        function set.IQsamples(obj, iSamples)
            if isa(iSamples,obj.getSavingFormat()) && size(iSamples,2) == 2
                obj.IQsamples = iSamples;
            else
                disp("Error, trying to return a wrong IQsamples array")
            end
        end

        %                        %
        % PRIVATE METHODS        %
        %                        %
     
        function tmpData = openReadCloseBinaryFile(obj,currentSample,nSamples)
            try
                fid = fopen(obj.inFullfilename);
                if(fid~=-1)
                    fseek(fid, currentSample*obj.nByte_per_sample*2, 'bof');
                    tmpData = fread(fid,[2 nSamples],obj.getSampleFormat());
                else 
                    disp("Matlab cannot open the file: aborted",obj.inFullfilename);
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

        function format = getSavingFormat(obj)
            if obj.nByte_per_sample == 1
                format = "int8";
            elseif obj.nByte_per_sample == 2
                format = "int16";
            elseif obj.nByte_per_sample == 3
                format = "bit24";
            elseif obj.nByte_per_sample == 4
                format = "int32";
            else
                disp("Errore: impossible value for nByte_per_sample")
            end
        end
        
   end
end