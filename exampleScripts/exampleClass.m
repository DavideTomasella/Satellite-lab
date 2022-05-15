%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% CLASS USAGE %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INSTANCE AND CONFIGURATION
%------------------------------------
inInterface = InOutInterface();
dirName = "./inData";
inInterface.configReadBinaryFile(dirName);

% PROCESS
%------------------------------------
filename = "data1.bin";
rawData = inInterface.readBinaryFile(filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef InOutInterface < handle
    %InOutInterface handles the acquisition of configuration and sampled
    %data and the output of the binary message
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
        %isGetConfigConfigured: transform the input rawData into a config
        %struct
        isGetConfigConfigured=false
        %isGetDataArrayConfigured: transform the input rawData into a
        %sampled data array
        isGetDataArrayConfigured=false
        %isReadBinaryFileConfigured: open a binary file and get its rawData
        isReadBinaryFileConfigured=false
        %isReadTxtFileConfigured: open a binary file and get its rawData
        isReadTxtFileConfigured=false
        %isSaveDataArrayConfigured: tranform the output data array into a
        %binary stream
        isSaveDataArrayConfigured=false
        %isWriteBinaryFileConfigured: write a binary stream into a file
        isWriteBinaryFileConfigured=false
    end
    properties (Constant)
        v="1"
    end
    properties(Dependent=true)
    end
    properties
       myDirectory        
    end

    methods
        function obj = InOutInterface()
            %InOutInterface constructor
            %   Create InOutInterface class
        end

        %                        %
        % CONFIGURATION METHODS  %
        %                        %

        function obj = configReadBinaryFile(obj,dirName)
            %configReadBinaryFile config file location, how to open the files,
            %and how to acquire the data
            %   dirName: string with absolute or relative path where
            %   configuration file and sampled data are stored
            if nargin < 2
                dirName = "./";
            end
            obj.myDirectory = dir(dirName);
            obj.isReadBinaryFileConfigured = true;
        end

        %                        %
        % ACTION METHODS         %
        %                        %

        function rawData = readBinaryFile(obj,filename)
            %readBinaryFile open the given file and return the bitstream of
            %its content
            %   filename: string with the name of the input file
            rawData=0;
            if obj.isReadBinaryFileConfigured
                rawData=0;
            end
        end

        %                        %
        % PRIVATE METHODS        %
        %                        %

    end
end