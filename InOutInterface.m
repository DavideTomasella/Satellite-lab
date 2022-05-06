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
        %isGetSettingsConfigured: transform the json input file into a config
        %struct and validate the input parameters
        isGetSettingsConfigured=false
        %isSaveResultsConfigured: write a file from a structure with the
        %demodulation results
        isSaveResultsConfigured=false
    end
    properties (Constant)
        v="1"
    end
    properties(Dependent=true)
    end
    properties
       inDirectory
       outDirectory        
    end

    methods
        function obj = InOutInterface()
            %InOutInterface constructor
            %   Create InOutInterface class
        end

        %                        %
        % CONFIGURATION METHODS  %
        %                        %

        function obj = configGetSettings(obj,dirName)
            %configGetSettings transform the json input file into a config
            %struct and validate the input parameters
            %   dirName: string with absolute or relative path where
            %   configuration file and sampled data are stored
            if nargin < 2
                dirName = "./";
            end
            obj.inDirectory = dir(dirName);
            obj.isGetSettingsConfigured = true;
        end

        function obj = configSaveResults(obj,dirName)
            %configSaveResults config write a file from a structure with
            %the demodulation results
            %   dirName: string with absolute or relative path where
            %   configuration file and sampled data are stored
            if nargin < 2
                dirName = "./";
            end
            obj.outDirectory = dir(dirName);
            obj.isGetSettingsConfigured = true;
        end
        
        %                        %
        % ACTION METHODS         %
        %                        %

        function obj = getSettings(obj)
            %getSettings transform the json input file into a config
            %struct and validate the input parameters
            %   ...
        end

        function obj = saveResults(obj)
            %saveResults write a file from a structure with the
            %demodulation results
            %   ...
        end
        
        %                        %
        % PRIVATE METHODS        %
        %                        %

        function obj = readJsonFile(obj)
            %readJsonFile read the content of a json file into a struct
            %   ...
        end

        function obj = validateSettings(obj)
            %validateSettings validate the input setting parameters
            %   ...
        end

        function obj = writeJsonFile(obj)
            %writeJsonFile write a json file from a structure
            %   ...
        end

    end
end