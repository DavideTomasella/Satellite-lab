%
% First implementation: Lorenzo Borsoi & Davide Tomasella (Tracking algorithm)
%                       Lorenzo Borsoi (Message Analyzer)
% Review and Testing: Davide Tomasella & Lorenzo Borsoi
%
classdef TT_Demodulator < handle %%TODO TrackingDemodulator
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

        function [decSymbols, idShift, idDoppler] = decodeOptimumShift(obj, inSamples, segmentSize, ...
                                                                       shifts_delayPRN, shifts_nSamples_x_symbolPeriod, ...
                                                                       shifts_nSamples_x_chipPeriod, coherenceFraction)    
            %windowSize: numero di simboli
            %nchips: numero di chips che compongono la PNR sequence
            %numero di PNR sequence in un simbolo: symbolPeriod/(chipPeriod*nchips)
            %numero di PNR sequence in M messaggi: windowSize*symbolPeriod/(chipPeriod*nchips)
            %PNRsequence va upsampled con fsampling in PNRcodes 1xN

            %Reshape input to have predetermined vector
            if size(inSamples, 1) < size(inSamples, 2)
                inSamples = inSamples'; %column vector
            end
            if size(shifts_nSamples_x_chipPeriod, 1) < size(shifts_nSamples_x_chipPeriod, 2)
                shifts_nSamples_x_chipPeriod = shifts_nSamples_x_chipPeriod'; %column vector
            end
            if size(shifts_nSamples_x_symbolPeriod, 1) < size(shifts_nSamples_x_symbolPeriod, 2)
                shifts_nSamples_x_symbolPeriod = shifts_nSamples_x_symbolPeriod'; %column vector
            end
            if size(shifts_delayPRN, 2) < size(shifts_delayPRN, 1)
                shifts_delayPRN = shifts_delayPRN'; %row vector
            end
            %create PRN with different length due to different dopplers
            PRNsampled = obj.getPRNFromChipPeriods(shifts_nSamples_x_chipPeriod, segmentSize);

            %adapt the length of the samples to the PRN, horizontal vector
            mySamples = obj.adaptSamplesToPRNlength(inSamples, size(PRNsampled, 2));
            
            %add different delays to the PRN, sampled PRN with shifts in time and frequency
            shiftedPRNsampled = obj.createShiftedPRN(PRNsampled, shifts_delayPRN);
            %plot(shiftedPRNsampled(1:2,end-1e4:1:end)')
            
            %in-phase & quadrature multicorrelation
            corrI = obj.normMultiply(shiftedPRNsampled, mySamples(1, :));
            corrQ = obj.normMultiply(shiftedPRNsampled, mySamples(2, :));
            
            %coherent integration, over the rows there are the coherent sums
            coherentCorrI = obj.sumOverCoherentFraction(corrI, shifts_nSamples_x_symbolPeriod, ...
                                                        coherenceFraction, segmentSize);
            coherentCorrQ = obj.sumOverCoherentFraction(corrQ, shifts_nSamples_x_symbolPeriod, ...
                                                        coherenceFraction, segmentSize);

            %non-coherent integration, column vector
            noncoherentCorr = sum(coherentCorrI .^ 2 + coherentCorrQ .^ 2, 2);

            %find max
            [~, idMax] = max(noncoherentCorr,[], 1);           
            [idDoppler, idShift] = ind2sub([size(shifts_nSamples_x_chipPeriod,1) size(shifts_delayPRN,2)], idMax);         
            %NOTE: DOT product: sommatoria[(Ie-Il)*Ip] - sommatoria[(Qe-Ql)*qp],
            %if <0 detector is late, if >0 too early

            %complete correlation over symbols
            %TODO channel inversion? linear estimator?
            bestCorrI = obj.sumFractionsOverSymbols(coherentCorrI(idMax, :), coherenceFraction);
            bestCorrQ = obj.sumFractionsOverSymbols(coherentCorrQ(idMax, :), coherenceFraction);
            
            %decoding
            decSymbols = (2 * (bestCorrI > 0) - 1)'; % column vector of decoded symbols +1,-1            
            
        end

        function PRNsampled = getPRNFromChipPeriods(obj, shifts_nSamples_x_chipPeriod, windowSize)
            optionDopplerLength = size(shifts_nSamples_x_chipPeriod, 1);
            maxLength = max(shifts_nSamples_x_chipPeriod * length(obj.PRNsequence) * ...
                            obj.nPRN_x_Symbol * windowSize);
            PRNsampled = zeros(optionDopplerLength, int32(maxLength + 0.5)); %ceil approx

            for period = 1:optionDopplerLength
                PRNinterp = interp1(1:length(obj.PRNsequence), obj.PRNsequence, ...
                    1:1 / shifts_nSamples_x_chipPeriod(period):length(obj.PRNsequence), "previous"); %upsampling   
                repSize = size(PRNinterp,2) * windowSize * obj.nPRN_x_Symbol;
                PRNsampled(period, 1:repSize) = repmat(PRNinterp, 1, windowSize * obj.nPRN_x_Symbol);
            end
        end

        function  mySamples = adaptSamplesToPRNlength(~, filteredSamples, PRNlength)
            if size(filteredSamples, 1) < PRNlength
                mySamples = [filteredSamples; zeros(PRNlength - size(filteredSamples, 1), ...
                                                    PRNlength)]';
            else
                mySamples = filteredSamples(1:PRNlength, :)';
            end
        end

        function shiftedPRNsampled = createShiftedPRN(~, PRNsampled, shifts_delayPRN)
            optionDopplerLength = size(PRNsampled, 1);
            optionDelayLength = size(shifts_delayPRN, 2);
            shiftedPRNsampled = zeros(optionDelayLength * optionDopplerLength, ...
                                      size(PRNsampled, 2));
            for shiftID = 0:optionDelayLength - 1
                offset = shiftID * optionDopplerLength;
                shiftedPRNsampled(1 + offset:offset + optionDopplerLength, :) = ...
                            circshift(PRNsampled, shifts_delayPRN(1 + shiftID), 2);
            end
        end

        function corr = normMultiply(~, shiftedPRNsampled, mySamples)
            corr = shiftedPRNsampled .* mySamples ./ sum(abs(shiftedPRNsampled) > 0, 2);
        end

        function coherentCorr = sumOverCoherentFraction(~, corr, shifts_nSamples_x_symbolPeriod, ...
                                                         coherenceFraction, windowSize)
            nCoherentSegments = windowSize / coherenceFraction;
            optionLength = size(corr, 1);
            corrLength = size(corr, 2);
            %row vector of coherent sum, option in each row
            coherentCorr = zeros(optionLength, nCoherentSegments);
            for cLine = 0:optionLength - 1 %cycle on columns
                symPeriod = shifts_nSamples_x_symbolPeriod(1 + mod(cLine, size(shifts_nSamples_x_symbolPeriod, 1)));
                %DT$ TODO check this approximation
                cohSamp = int32(symPeriod * coherenceFraction);
                cLine_padded = zeros(int32(corrLength / cohSamp + 0.5) * cohSamp,1);
                cLine_padded(1:corrLength) = corr(1 + cLine,:)';
                corrSum = sum(reshape(cLine_padded, cohSamp, []), 1);
                coherentCorr(1 + cLine,:) = corrSum(1:nCoherentSegments);
            end
        end

        function bestCorr = sumFractionsOverSymbols(~, bestCoherentCorr, coherenceFraction)
            %columns of coherent fractions for each symbol
            tmpCorr = reshape(bestCoherentCorr, 1 / coherenceFraction, []);
            bestCorr = sum(tmpCorr, 1);
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