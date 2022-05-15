%
% First implementation: Lorenzo Borsoi
% Review and Testing: Davide Tomasella
%
classdef Demodulator < handle
    %Demodulator handles...

    properties 

        fsampling
        PRNsequence
        nPRN_x_Symbol
        CRCkey
        PRN_SV_ID
        SyncPattern
        SV_IDLength
        M_IDLength 
        M_bodyLength
        CRCLength 
        TailPattern
        %windowSize
        %symbolPeriod
        %chipPeriod
        %coherenceFraction
    end

    methods
        function obj = Demodulator()
            %Demodulator constructor
            %   Create Demodulator class
        end

        function obj = configCorrelatorValues(obj,fSampling, nPRN_x_Symbol, PRNcode)
            obj.fsampling = fSampling;
            obj.nPRN_x_Symbol = nPRN_x_Symbol;
            obj.PRNsequence = 2 * PRNcode - 1;   
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

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  SYMBOLS TRACKING                  $
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [decSymbols, idShift, idDoppler] = decodeOptimumShift(obj, filteredSamples, windowSize, ...
                                                                       sampleShifts, e_nSamples_x_symbolPeriod, ...
                                                                       e_nSamples_x_chipPeriod, coherenceFraction)    
            %windowSize: numero di simboli
            %wsize: numero campioni
            %nchips: numero di chips che compongono la PNR sequence
            %numero di PNR sequence in un simbolo: symbolPeriod/(chipPeriod*nchips)
            %numero di PNR sequence in M messaggi: windowSize*symbolPeriod/(chipPeriod*nchips)
            %PNRsequence va upsampled con fsampling in PNRcodes 1xN

            %Reshape input to have predetermined vector
            if size(e_nSamples_x_chipPeriod, 1) < size(e_nSamples_x_chipPeriod, 2)
                e_nSamples_x_chipPeriod = e_nSamples_x_chipPeriod'; %column vector
            end
            if size(e_nSamples_x_symbolPeriod, 1) < size(e_nSamples_x_symbolPeriod, 2)
                e_nSamples_x_symbolPeriod = e_nSamples_x_symbolPeriod'; %column vector
            end
            if size(sampleShifts, 2) < size(sampleShifts, 1)
                sampleShifts = sampleShifts'; %row vector
            end
            PRNsampled = obj.getPRNsampledFromWindow(windowSize, e_nSamples_x_chipPeriod);

            %mySamples horizontal vector
            if size(filteredSamples, 1) < size(PRNsampled, 2)
                mySamples = [filteredSamples; zeros(size(PRNsampled, 2) - size(filteredSamples, 1), ...
                                                    size(filteredSamples, 2))]';
            else
                mySamples = filteredSamples(size(PRNsampled, 2), :)';
            end
            
            %sampled PRN with given shifts in time and frequency            
            shiftedPRNsampled = zeros(size(sampleShifts, 2) * size(e_nSamples_x_chipPeriod, 1), ...
                                      size(PRNsampled, 2));
            for shiftID = 0:size(sampleShifts, 2) - 1
                offset = shiftID * size(e_nSamples_x_chipPeriod, 1);
                shiftedPRNsampled(1 + offset:offset + size(e_nSamples_x_chipPeriod, 1), :) = ...
                            circshift(PRNsampled, sampleShifts(1 + shiftID), 2);
            end
            %plot(shiftedPRNsampled(1:2,end-1e4:1:end)')
            
            %in-phase multicorrelation
            corrI = shiftedPRNsampled .* mySamples(1, :) ./ sum(abs(shiftedPRNsampled) > 0, 2);
            %quadrature multicorrelation
            corrQ = shiftedPRNsampled .* mySamples(2, :) ./ sum(abs(shiftedPRNsampled) > 0, 2);
            
            %coherent integration
            nCoherentSegments = windowSize / coherenceFraction;
            coherentCorrI = zeros(size(corrI, 1), nCoherentSegments);%row vector of coherent sum
            coherentCorrQ = zeros(size(corrQ, 1), nCoherentSegments);
            corrLength = size(corrI, 2);
            for cI = 0:size(corrI, 1) - 1 %cycle on column
                symPeriod = e_nSamples_x_symbolPeriod(1 + fix(cI / size(sampleShifts, 2)));
                cohSamp = int32(symPeriod * coherenceFraction);
                cI_padded = zeros(int32(corrLength / cohSamp + 0.5) * cohSamp,1);
                cI_padded(1:corrLength) = corrI(1 + cI,:)';
                corrSum = sum(reshape(cI_padded, cohSamp, []), 1);
                coherentCorrI(1 + cI,:) = corrSum(1:nCoherentSegments);
            end
            for cQ = 0:size(corrQ, 1) - 1 %cycle on column
                symPeriod = e_nSamples_x_symbolPeriod(1 + fix(cQ / size(sampleShifts, 2)));
                cohSamp = int32(symPeriod * coherenceFraction);
                cQ_padded = zeros(int32(corrLength / cohSamp + 1) * cohSamp,1);
                cQ_padded(1:corrLength) = corrQ(1 + cQ,:)';
                corrSum = sum(reshape(cQ_padded, cohSamp, []), 1);
                coherentCorrQ(1 + cQ,:) = corrSum(1:nCoherentSegments);
            end

            %non-coherent integration column vector
            noncoherentCorr = sum(coherentCorrI .^ 2 + coherentCorrQ .^ 2, 2);

            %find max
            [~, idMax] = max(noncoherentCorr,[], 1);           
            [idDoppler, idShift] = ind2sub([size(e_nSamples_x_chipPeriod,1) size(sampleShifts,2)], idMax);         
            %NOTE: DOT product: sommatoria[(Ie-Il)*Ip] - sommatoria[(Qe-Ql)*qp],
            %if <0 detector is late, if >0 too early

            %decoding
            bestCorrI = reshape(coherentCorrI(idMax, :), 1 / coherenceFraction, []); %columns of coherent fractions for each symbol
            bestCorrQ = reshape(coherentCorrQ(idMax, :), 1 / coherenceFraction, []); %columns of coherent fractions for each symbol
            decSymbols = (2 * (sum(bestCorrI, 1) > 0) - 1)'; % column vector of decoded symbols +1,-1            
            
        end

        function PRNsampled = getPRNsampledFromWindow(obj, windowSize, e_nSamples_x_chipPeriod)
            maxLength = max(e_nSamples_x_chipPeriod * length(obj.PRNsequence) * ...
                            obj.nPRN_x_Symbol * windowSize);
            PRNsampled = zeros(size(e_nSamples_x_chipPeriod,1), int32(maxLength + 0.5)); %ceil approx

            for period = 1:size(e_nSamples_x_chipPeriod,1)
                PRNinterp = interp1(1:length(obj.PRNsequence), obj.PRNsequence, ...
                    1:1 / e_nSamples_x_chipPeriod(period):length(obj.PRNsequence), "previous"); %upsampling   
                repSize = size(PRNinterp,2) * windowSize * obj.nPRN_x_Symbol;
                PRNsampled(period, 1:repSize) = repmat(PRNinterp, 1, windowSize * obj.nPRN_x_Symbol);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  MESSAGE DECODING                  $
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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