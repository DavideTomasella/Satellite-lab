%
% First implementation: Davide Tomasella
% Review and Testing:
%
classdef InOutInterface < handle
    %InOutInterface handles the acquisition of configuration and the
    %saving of the output file (demodulation results)
    %   The class usage is performed through 3 steps:
    %   1- class initialization (without particular arguments)
    %   2- class configuration through the configACTION methods that get in
    %   input the configuration parameters and defines the ACTION procedure
    %   from them. These methods must not save array of data as internal
    %   variables to minimize the space usage. At the end of the
    %   configuration procedure a protected state isACTIONConfigured=true
    %   must be set.
    %   3- class operation through ACTION methods that performed on the
    %   input data array the ACTION previously configured in configACTION
    %   methods and the value of local variables. The ACTION algorithm must
    %   be performed only after isACTIONConfigured is set.

    properties (Access = { ?classUnderTest, ?matlab.unittest.TestCase })
        %isCreateSettingsConfigured: transform the json input file into a config
        %struct and validate the input parameters
        isCreateSettingsConfigured=false
        %isSaveResultsConfigured: write a file from a structure with the
        %demodulation results
        isSaveResultsConfigured=false
    end
    properties (Constant)
        v="1"
    end
    properties(SetAccess=immutable, GetAccess=public)
       defaultSettings
    end
    properties (SetAccess=private, GetAccess=public)
       inDirectory
       outDirectory
       settings
       PRNcode
       settings_lastLoad_filename
       PRNcode_lastLoad_filename
       results_lastSave_filename
    end
    properties (SetAccess=public, GetAccess=public)
        results
    end

    methods
        function obj = InOutInterface()
            %InOutInterface constructor
            %   Create InOutInterface class
            
            obj.defaultSettings = obj.createDefaultSettings();
            obj.results = obj.createResultsPlaceholder();
        end

        %                        %
        % CONFIGURATION METHODS  %
        %                        %

        function obj = configCreateSettings(obj,dirName)
            %configCreateSettings transform the json input file into a config
            %struct and validate the input parameters
            %   dirName: string with absolute or relative path where
            %   configuration file and sampled data are stored

            if nargin < 2
                dirName = "./";
            end
            obj.inDirectory = dir(dirName);
            obj.isCreateSettingsConfigured = true;
        end

        function obj = configSaveResults(obj,dirName,relativePath)
            %configSaveResults config write a file from a structure with
            %the demodulation results
            %   dirName: string with absolute or relative path where
            %   configuration file and sampled data are stored
            %   relativePath: true if the dirName is relative to the
            %   current directory

            if nargin < 3
                relativePath = true;
            end
            if nargin < 2
                dirName = "./";
            end

            if relativePath
                dirName = fullfile(".",dirName);
            end
            obj.outDirectory = dir(dirName);
            obj.isCreateSettingsConfigured = true;
        end
        
        %                        %
        % ACTION METHODS         %
        %                        %

        function oSettings = createSettings(obj,filename,prnfilename)
            %createSettings transform the json input file into a config
            %struct and validate the input parameters
            %   filename:name of input json file
            %   prnfilename:name of input json file with prn sequences

            if nargin>=2 && obj.isCreateSettingsConfigured
                try %Read settings file
                    [newSettings, newFilename] = obj.readJsonFile(filename);
                    try %Validate settings
                        if obj.validateSettings(newSettings)
                            if nargin>=3
                                %Save new settings
                                obj.settings = newSettings;
                                obj.settings_lastLoad_filename = newFilename;
                                try %Read all PRN and select one
                                    [newPRNlist, newPrnfilename] = obj.readJsonFile(prnfilename);
                                    %disp(newPRNlist.("E1B_Code_No_"+newSettings.SV_PRN_ID))
                                    newPRNcode = hexToBinaryVector(newPRNlist.("E1B_Code_No_"+newSettings.SV_PRN_ID));
                                catch
                                    disp("PRN pattern not valid: TO HANDLE", newSettings.PRNcode);
                                end
                                %Save new PRN
                                obj.PRNcode = newPRNcode;
                                obj.PRNcode_lastLoad_filename = newPrnfilename;
                            end
                        else
                            disp("new settings not valid: use default", newSettings);
                        end                            
                    catch
                        disp("error in settings validation: use default", newSettings,newFilename);
                    end
                catch
                    disp("error in json acquisition: use default");
                end
            end
            % Return the created settings or the default ones
            oSettings = obj.settings;
        end

        function oSettings = get.settings(obj)
            %getSettings returns the loaded settings of the default ones

            if ~isempty(obj.settings)
                if length(fieldnames(obj.settings)) >= length(fieldnames(obj.defaultSettings))
                    oSettings = obj.settings;    
                    return
                end
            end
            oSettings = obj.defaultSettings;
        end

        function oFilename = saveResults(obj,filename)
            %saveResults write a file from a structure with the
            %demodulation results adding date/hour to the name
            %   filename: output json filename

            if nargin < 2
                disp("specify json filename: aborted");
            end
            try
                [subdir,name,ext] = fileparts(filename);
                date = datestr(now,'_yymmdd_HHMMSS');
                newFilename = strcat(subdir,name,date,ext);
                obj.writeJsonFile(newFilename);
                obj.results_lastSave_filename = newFilename;
            catch
                disp("impossible to save json file: aborted", newFilename);
            end
            % Return the created filename or null
            oFilename = obj.results_lastSave_filename;
        end
        
        %                        %
        % PRIVATE METHODS        %
        %                        %
        
        function defSettings = createDefaultSettings(~)
            %createDefaultSettings default settings struct

            defSettings = struct("fSampling", 10e6,...
                                 "quantizationBits", 16,...
                                 "scenarioDuration", 5.05,...
                                 "SV_PRN_ID", 1,...
                                 "CRCpolynomial", "A23DCB",...
                                 "SYNCpattern", "0101100000",...
                                 "TAILpattern", "000000",...
                                 "SVIDlength", 6,...
                                 "MIDlength", 4,...
                                 "MBODYlength_TX", 80,...
                                 "MBODYlength_ACK", 30,...
                                 "CRClength", 24,...
                                 "nPRN_x_Symbol", 1,...
                                 "nChip_x_PRN", 4092,...
                                 "chipRate", 1.023e6);
        end

        function dddResults = createResultsPlaceholder(~)
            %createResultsPlaceholder results struct placeholder

            dddResults = struct("version","1", ...
                                "SV_ID","000000",...
                                "message_ID","0000",...
                                "message_body","000000000000000000000000000000",...
                                "CRC","000000000000000000000000",...
                                "ACKed",false,...
                                "isACKmessage",false,...
                                "estimatedDoppler",0.0,...
                                "estimatedDelay",0.0);
        end

        function [jsonStruct,fullname] = readJsonFile(obj,filename)
            %readJsonFile read the content of a json file into a struct
            %   filename:name of input json file

            folder = obj.inDirectory.folder;
            fullname = fullfile(folder,filename);
            jsonStruct = jsondecode(fileread(fullname));
        end

        function isValid = validateSettings(obj,newSettings)
            %validateSettings validate the input setting parameters
            %comparing it with default one
            %   newSettings: setting struct to validate
            try
                for setName = fieldnames(obj.defaultSettings)'
                    newSettings.(string(setName));
                end
                isValid = true;
            catch
                disp("Error: new settings lacks at least one parameter")
            isValid = false;
            end
        end

        function jsonText = writeJsonFile(obj,filename)
            %writeJsonFile write a json file from a structure
            %   filename:name of output json file

            jsonText = jsonencode(obj.results, "PrettyPrint", true);
            folder = obj.outDirectory.folder;
            fullname = fullfile(folder,filename);
            fid=fopen(fullname,'wt');
            fprintf(fid, jsonText);
            fclose(fid);
        end

    end
end