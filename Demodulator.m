%
% First implementation: Lorenzo Borsoi
% Review and Testing:
%
classdef Demodulator < handle
    %Demodulator handles...

    properties (SetAccess=private, GetAccess=public)

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
            
            prodToBin = @(a, b) dec2binvec(a * hex2dec(b));
            % (X + 1) * P(X)
            obj.CRCkey = prodToBin(3, CRCpolynomial); %test with actual CRC and messages  
            % same as [0 CRCpolynomial]+[CRCpolynomial 0]
            % TEST
            % if (binvec2dec(prodToBin(3,CRCpolynomial))-binvec2dec(hexToBinaryVector(CRCpolynomial)))/2-...
            %    binvec2dec(hexToBinaryVector(CRCpolynomial)) ~= 0 disp("Error"); end

            obj.PRN_SV_ID = PRN_SV_ID;
            obj.SyncPattern = SYNCpattern;
            obj.SV_IDLength = SVIDlength;
            obj.M_IDLength = MIDlength; 
            obj.M_bodyLength = MBODYlength; %verify if first bit is zero
            obj.CRCLength = CRClength; %compute and verify
            obj.TailPattern = TAILpattern; %%000000
            
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
            shiftedPRNsampled = zeros(size(sampleShifts,2) * size(e_nSamples_x_chipPeriod,1), ...
                                      size(PRNsampled,2));
            for shiftID = 1:size(sampleShifts,2)
                offset = (shiftID - 1) * size(e_nSamples_x_chipPeriod,1);
                shiftedPRNsampled(1 + offset:offset + size(e_nSamples_x_chipPeriod,1), :) = ...
                            circshift(PRNsampled, sampleShifts(shiftID), 2);
            end
            
            %in-phase multicorrelation
            corrI = shiftedPRNsampled .* mySamples(:, 1) / sum(abs(shiftedPRNsampled) > 0, 2);
            %quadrature multicorrelation
            corrQ = shiftedPRNsampled .* mySamples(:, 2) / sum(abs(shiftedPRNsampled) > 0, 2);
            
            %coherent integration
            coherentCorrI = zeros(1,size(corrI, 1));%row vector of coherent sum
            coherentCorrQ = zeros(1,size(corrQ, 1));
            for cI = 1:size(corrI,1)%cycle on column
                %TODO check padding
                symPeriod = e_nSamples_x_symbolPeriod(mod(cI, size(sampleShifts, 2)));
                cohSamp = int32(symPeriod * coherenceFraction);
                cI_padded = [corrI(cI,:)'; zeros(int32(length(cI) / cohSamp + 1) * cohSamp - length(cI), 1)];
                coherentCorrI(cI,:) = sum(reshape(cI_padded, cohSamp, []), 1);
            end
            for cQ = 1:size(corrQ,1)%cycle on column
                cohSamp = e_nSamples_x_symbolPeriod(mod(cQ, size(sampleShifts, 2)));
                cQ_padded = [corrQ(cQ,:)'; zeros(int32(length(cQ) / cohSamp + 1) * cohSamp - length(cQ), 1)];
                coherentCorrQ(cQ,:) = sum(reshape(cQ_padded, cohSamp, []), 1);
            end

            %non-coherent integration column vector
            noncoherentCorr = sum(coherentCorrI .^ 2 + coherentCorrQ .^ 2, 2);

            %find max
            [~, idMax] = max(noncoherentCorr, 1);           
            [idDoppler, idShift] = ind2sub([size(sampleShifts,2) size(e_nSamples_x_chipPeriod,1)], idMax);

            %decoding
            bestCorrI = reshape(coherentCorrI(idMax, :), 1 / coherenceFraction, []); %columns of coherent fractions for each symbol
            bestCorrQ = reshape(coherentCorrQ(idMax, :), 1 / coherenceFraction, []); %columns of coherent fractions for each symbol
            decSymbols = (2 * (sum(bestCorrI, 1) > 0) - 1)'; % column vector of decoded symbols +1,-1            
            
        end

        function PRNsampled = getPRNsampledFromWindow(windowSize, e_nSamples_x_chipPeriod)
            maxLength = max(e_nSamples_x_chipPeriod * length(obj.PRNsequence) * ...
                            obj.nPRN_x_Symbol * windowSize);
            PRNsampled = zeros(size(e_nSamples_x_chipPeriod,1), maxLength);
            for period = 1:size(e_nSamples_x_chipPeriod,1)
                PRNinterp = interp1(1:length(obj.PRNsequence), obj.PRNsequence, ...
                    1:1 / e_nSamples_x_chipPeriod:length(obj.PRNsequence), "previous"); %upsampling            
                PRNsampled(period, :) = repmat(PRNinterp, 1, windowSize * obj.nPRN_x_Symbol);
            end
        end

            %somma interna ad un coherence time -> coherent integration
            %somma dei quadrati delle somme risultanti -> non-coherent integration
            %discriminator
            %DOT product: sommatoria[(Ie-Il)*Ip] - sommatoria[(Qe-Ql)*qp],
            %if <0 detector is late, if >0 too early

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  MESSAGE DECODING                  $
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function analyzeMessage(obj, decodedSymbols, outData)

            %minimum distance region decision, no channel inversion
            decodedSymbols = int16((decodedSymbols + 1) ./ 2); % 1 -> 1   -1 -> 0

            [SV_ID, M_ID, M_body, CRCMessage] = obj.validateAndSplitMessage(decodedSymbols);

            outData.SV_ID = num2str(SV_ID);
            outData.message_ID = num2str(M_ID);    
            outData.message_Body = num2str(M_body);         
            outData.CRC = num2str(CRCMessage);                      
            
            CRCCheck = obj.calcChecksum([M_ID M_body CRCMessage]);   

            %verifications
            outData.isACKmessage = (M_body(1) == 0); 
            outData.ACKed = (sum(CRCCheck) == 0);
            if ~outData.ACKed
                sprintf("Error no ack: received %s, calculated %s", string(CRCmessage), string(CRCCheck))
            end
        end

        function [SV_ID, M_ID, M_body, CRCMessage] = validateAndSplitMessage(obj, decodedSymbols)
            messageLength = length(obj.SyncPattern) + obj.SV_IDLength + obj.M_IDLength + ...
                            obj.M_bodyLength + obj.CRCLength + length(obj.TailPattern);
            if length(decodedSymbols) < messageLength
                disp("error")
            end
            startPoint = 1;
            [SyncMessage, startPoint] = extractAndAdvance(decodedSymbols, startPoint, length(obj.SyncPattern));
            [SV_ID, startPoint] = extractAndAdvance(decodedSymbols, startPoint, obj.SV_IDLength);
            [M_ID, startPoint] = extractAndAdvance(decodedSymbols, startPoint, obj.M_IDLength);
            [M_body, startPoint] = extractAndAdvance(decodedSymbols, startPoint, obj.M_bodyLength);
            [CRCMessage, startPoint] = extractAndAdvance(decodedSymbols, startPoint, obj.CRCLength);
            [TailMessage, ~] = extractAndAdvance(decodedSymbols, startPoint, length(obj.TailPattern));            

            if TailMessage ~= obj.TailPattern 
                disp("error tail"); 
            end

            if SyncMessage ~= obj.SyncPattern 
                disp("error sync"); 
            end
            
            if binvec2dec(SV_ID) ~= obj.PRN_SV_ID
                disp("not my message")
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

            %padding
            tmpMessage = [messageToCheck, zeros(1, obj.CRCLength)]; %M_ID + M_body + CRC + padding

            j = find(tmpMessage, 1); %find symbol 1
            while sum(tmpMessage) > 0 && 1 + length(tmpMessage) - j > obj.CRCLength                   
                tmpMessage = tmpMessage(j:end);
                tmpMessage(1:1 + obj.CRCLength) = bitxor(tmpMessage(1:1 + obj.CRCLength), obj.CRCkey);                
                j = find(tmpMessage, 1);                
            end

            CRCCheck = tmpMessage(1 + end - obj.CRCLength: end);
        end

        function [newArray, nextPoint] = extractAndAdvance(~, originalArray, startPoint, length)
            %prende l'intervallo startPoint-startPoint+length e aggiorna startPoint+=length
            nextPoint = startPoint + length;
            newArray = originalArray(startPoint:nextPoint - 1);
        end
                    
    end
end