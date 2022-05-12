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
                                                                       sampleShifts, nSamples_x_symbolPeriod, ...
                                                                       e_nSamples_x_chipPeriod, coherenceFraction)    
            %windowSize: numero di simboli
            %wsize: numero campioni
            %nchips: numero di chips che compongono la PNR sequence
            %numero di PNR sequence in un simbolo: symbolPeriod/(chipPeriod*nchips)
            %numero di PNR sequence in M messaggi: windowSize*symbolPeriod/(chipPeriod*nchips)
            %PNRsequence va upsampled con fsampling in PNRcodes 1xN
            if size(e_nSamples_x_chipPeriod, 1) < size(e_nSamples_x_chipPeriod, 2)
                e_nSamples_x_chipPeriod = e_nSamples_x_chipPeriod'; %column vector
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
            shiftedPRNsampled = zeros(length(sampleShifts) * length(e_nSamples_x_chipPeriod), ...
                                      size(PRNsampled,2));
            for shift = sampleShifts
                shiftedPRNsampled(shift, :) = circshift(PRNsampled, shift, 2);
            end
            
            %in-phase multicorrelation
            corrI = shiftedPRNsampled .* mySamples(:, 1);
            %quadrature multicorrelation
            corrQ = shiftedPRNsampled .* mySamples(:, 2);
            
            %coherent integration
            coherentCorrI = zeros(size(corrI, 1),1);
            coherentCorrQ = zeros(size(corrQ, 1),1);
            for cI=corrI'%cycle on column
                %TODO check padding
                symPeriod = nSamples_x_symbolPeriod(mod(cI, length(sampleShifts)));
                cohSamp = int32(symPeriod*coherenceFraction);
                cI_padded = [cI'; int32(length(cI) / cohSamp + 1) * cohSamp - length(cI)];
                coherentCorrI(cI) = sum(reshape(cI_padded, cohSamp, []), 1);
            end
            for cQ = corrQ'%cycle on column
                cohSamp = nSamples_x_symbolPeriod(mod(cQ, length(sampleShifts)));
                cQ_padded = [cQ'; int32(length(cQ) / cohSamp + 1) * cohSamp - length(cQ)];
                coherentCorrQ(cQ) = sum(reshape(cQ_padded, cohSamp, []), 1);
            end

            %non-coherent integration
            noncoherentCorr = sum(coherentCorrI .^ 2 + coherentCorrQ .^ 2, 2);

            %find max
            [~, idMax] = max(noncoherentCorr, 1); % early -> 1 prompt -> 2 late -> 3            
            [idDoppler, idShift] = ind2sub([length(sampleShifts) length(e_nSamples_x_chipPeriod)], idMax);

            %decoding
            bestCorrI = reshape(coherentCorrI(idMax, :), 1 / coherenceFraction, []); %columns of coherent fractions for each symbol
            bestCorrQ = reshape(coherentCorrQ(idMax, :), 1 / coherenceFraction, []); %columns of coherent fractions for each symbol
            decSymbols = (2 * (sum(bestCorrI, 1) > 0) - 1)'; % column vector of decoded symbols +1,-1            
            
        end

        function PRNsampled = getPRNsampledFromWindow(windowSize, e_nSamples_x_chipPeriod)
            maxLength = max(e_nSamples_x_chipPeriod * length(obj.PRNsequence) * ...
                            obj.nPRN_x_Symbol * windowSize);
            PRNsampled = zeros(length(e_nSamples_x_chipPeriod), maxLength);
            for period = e_nSamples_x_chipPeriod
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
                disp("no ack", CRCmessage, CRCCheck)
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


        function CRCCheck = calcChecksum (obj, messageToCheck)
            %leftmost symbol equal to 1
            %compute checksum only on the superposed portion 
            %CRC
            %G(X) = (X+1)P(X)
            %P(X) = X23 + X17 + X13 + X12 + X11 + X9 + X8 + X7 + X5 + X3 + 1
            %     = 1000 0010 0011 1011 1010 1001
            %     = A23DCB  

            %padding
            tmpMessage = [messageToCheck, zeros(1, obj.CRCLength)]; %M_ID + M_body + CRC + padding

            j = find(tmpMessage, 1); %find symbol 1            
            while j > 0 && 1 + length(tmpMessage) - j > obj.CRCLength                   
                tmpMessage = tmpMessage(j:end);
                tmpMessage(1:1 + obj.CRCLength) = bitxor(tmpMessage(1:1 + obj.CRCLength), obj.CRCkey);                
                j = find(tmpMessage, 1);                
            end

            CRCCheck = tmpMessage(1 + end - obj.CRCLength: end);
        end

        function [newArray, startPoint] = extractAndAdvance(~, originalArray, startPoint, length)
            %prende l'intervallo startPoint-startPoint+length e aggiorna startPoint+=length

            newArray = originalArray(startPoint:startPoint + length - 1);             
            startPoint = startPoint+length;  
        end
                    
    end
end