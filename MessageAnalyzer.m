%
% First implementation: Lorenzo Borsoi
% Review and Testing: Lorenzo Borsoi
%
classdef MessageAnalyzer < handle
    %Demodulator handles...

    properties 
        CRCkey
        PRN_SV_ID
        SyncPattern
        SV_IDLength
        M_IDLength 
        M_bodyLength
        CRCLength 
        TailPattern
    end

    methods
        function obj = MessageAnalyzer()
            %MessageAnalyzer constructor
            %   Create MessageAnalyzer class
            
        end

        function obj = configMessageAnalyzer(obj, CRCpolynomial, PRN_SV_ID,...
                                             SYNCpattern, SVIDlength, MIDlength, ...
                                             MBODYlength, CRClength, TAILpattern)            
            
            % TEST
            % if (binvec2dec(prodToBin(3,CRCpolynomial))-binvec2dec(hexToBinaryVector(CRCpolynomial)))/2-...
            %    binvec2dec(hexToBinaryVector(CRCpolynomial)) ~= 0 disp("Error"); end

            try
                prodToBin = @(a, b) dec2binvec(a * hex2dec(b));
                % (X + 1) * P(X)
                obj.CRCkey = prodToBin(3, CRCpolynomial);
            catch
                log2dec = @(v) bin2dec(num2str(v));
                sumbin = @(a,b) dec2bin(sum([log2dec(a),log2dec(b)]))-'0';
                polArray = hexToBinaryVector(CRCpolynomial);
                obj.CRCkey = sumbin([0 polArray],[polArray 0]);
            end

            obj.PRN_SV_ID = PRN_SV_ID;
            obj.SyncPattern = SYNCpattern; %string "0101100000" to binary vector
            obj.SV_IDLength = SVIDlength;
            obj.M_IDLength = MIDlength; 
            obj.M_bodyLength = MBODYlength; %verify if first bit is zero
            obj.CRCLength = CRClength; %compute and verify
            obj.TailPattern = TAILpattern; % string "000000" to binary vector %char()-'0'
            
        end

        function analyzeMessage(obj, decodedSymbols, outInterface)

            if size(decodedSymbols, 1) > size(decodedSymbols, 2)
                decodedSymbols = decodedSymbols'; %column vector
            end            
            
            %minimum distance region decision, no channel inversion
            decodedSymbols = int16((decodedSymbols + 1) ./ 2); % 1 -> 1   -1 -> 0

            [SV_ID, M_ID, M_body, CRCMessage] = obj.validateAndSplitMessage(decodedSymbols);

            outInterface.results.SV_ID = num2str(SV_ID,'%d');
            outInterface.results.message_ID = num2str(M_ID,'%d');    
            outInterface.results.message_body = num2str(M_body,'%d');         
            outInterface.results.CRC = num2str(CRCMessage,'%d');             
            
            CRCCheck = obj.calcChecksum([M_ID M_body CRCMessage]);   

            %verifications
            outInterface.results.isACKmessage = M_body(1) == 0; 
            outInterface.results.ACKed = sum(CRCCheck) == 0;
            if ~outInterface.results.ACKed
                fprintf("Error no ack: received CRC %s and CRC remainder %s", string(num2str(CRCMessage,'%d')), string(num2str(CRCCheck,'%d')));
            end
        end

        function [SV_ID, M_ID, M_body, CRCMessage] = validateAndSplitMessage(obj, decodedSymbols)
            messageLength = length(obj.SyncPattern) + obj.SV_IDLength + obj.M_IDLength + ...
                            obj.M_bodyLength + obj.CRCLength + length(obj.TailPattern);
            if length(decodedSymbols) ~= messageLength
                warning('error! Message length is incorrect')
            end
            startPoint = 1;
            [SyncMessage, startPoint] = obj.extractAndAdvance(decodedSymbols, startPoint, length(obj.SyncPattern));
            [SV_ID, startPoint] = obj.extractAndAdvance(decodedSymbols, startPoint, obj.SV_IDLength);
            [M_ID, startPoint] = obj.extractAndAdvance(decodedSymbols, startPoint, obj.M_IDLength);
            [M_body, startPoint] = obj.extractAndAdvance(decodedSymbols, startPoint, obj.M_bodyLength);
            [CRCMessage, startPoint] = obj.extractAndAdvance(decodedSymbols, startPoint, obj.CRCLength);
            [TailMessage, ~] = obj.extractAndAdvance(decodedSymbols, startPoint, length(obj.TailPattern));            

            if string(num2str(TailMessage,'%d')) ~= string(obj.TailPattern) 
                warning('Error. Incorrect received tail %s.',num2str(TailMessage,'%d'))
            end

            if string(num2str(SyncMessage,'%d')) ~= string(obj.SyncPattern) 
                warning('Error. Incorrect received sync %s.',num2str(SyncMessage,'%d')); 
            end
            
            if bin2dec(char(SV_ID+'0')) ~= obj.PRN_SV_ID
                warning('Error. Not my message, incorrect received SV_ID %s.',char(SV_ID+'0'))
            end
        end

        function CRCCheck = calcChecksum(obj, messageToCheck)
            %leftmost symbol equal to 1
            %compute checksum only on the superposed portion 
            %CRC
            %G(X) = (X+1)P(X)
            %P(X) = X^23 + X^17 + X^13 + X^12 + X^11 + X^9 + X^8 + X^7 + X^5 + X^3 + 1
            %     = 1000 0010 0011 1011 1010 1001
            %     = A23DCB  

            %initialization
            tmpMessage = double(messageToCheck); %M_ID + M_body + CRC

            j = find(tmpMessage, 1); %find symbol 1
            while sum(tmpMessage) > 0 && 1 + length(tmpMessage) - j > obj.CRCLength                   
                tmpMessage = tmpMessage(j:end);
                tmpMessage(1:1 + obj.CRCLength) = bitxor(tmpMessage(1:1 + obj.CRCLength), obj.CRCkey);                
                j = find(tmpMessage, 1);                
            end

            CRCCheck = tmpMessage(1 + end - obj.CRCLength: end);
        end

        function [newArray, nextPoint] = extractAndAdvance(~, originalArray, startPoint, l)
            %prende l'intervallo startPoint-startPoint+length e aggiorna startPoint+=length
            nextPoint = startPoint + l;
            newArray = originalArray(startPoint:nextPoint - 1);
        end
                    
    end
end