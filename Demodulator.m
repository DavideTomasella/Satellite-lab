%
% First implementation: Lorenzo Borsoi
% Review and Testing:
%
classdef Demodulator < handle
    %Demodulator handles...

    properties (SetAccess=public, GetAccess=public)

        CRCkey
        SyncLength
        SV_IDLength
        M_IDLength 
        M_bodyLength
        CRCLength 
        TailLength         
        Sync
        Tail
        fsampling
        PRNsequence
        nPRN_x_Symbol

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

        function obj = configCorrelatorValues(obj,fSampling,nPRN_x_Symbol,PRNsequence)
            obj.fsampling=fSampling;
            obj.nPRN_x_Symbol=nPRN_x_Symbol;
            obj.PRNsequence = PRNsequence;   
        end

        function obj = configMessageAnalyzer(obj,CRCpolynomial, ...
                                            SYNCpattern,SVIDlength,MIDlength, ...
                                            MBODYlength,CRClength,TAILpattern)
                        
            log2dec = @(v) bin2dec(num2str(v));
            sumbin = @(a,b) dec2bin(sum([log2dec(a),log2dec(b)]))-'0';
            polArray = hexToBinaryVector(CRCpolynomial);
            % P(X) + X*P(X)            
            obj.CRCkey = sumbin([0 polArray],[polArray 0]); %test with actual CRC and messages        
            obj.Sync=SYNCpattern;
            obj.SyncLength=length(SYNCpattern);
            obj.SV_IDLength=SVIDlength;
            obj.M_IDLength=MIDlength; 
            obj.M_bodyLength=MBODYlength; %verify if first bit is zero
            obj.CRCLength=CRClength; %compute and verify
            obj.Tail=TAILpattern; %%000000
            obj.TailLength=length(TAILpattern); 
            
        end


        function [decSymbols,idShift] = decodeOptimumShift(obj,filteredSamples,windowSize, ...
                                                            sampleShifts,nSamples_x_symbolPeriod, ...
                                                            coherenceFraction)    
            %windowSize: numero di simboli
            %wsize: numero campioni
            %nchips: numero di chips che compongono la PNR sequence
            %numero di PNR sequence in un simbolo: symbolPeriod/(chipPeriod*nchips)
            %numero di PNR sequence in M messaggi: windowSize*symbolPeriod/(chipPeriod*nchips)
            %PNRsequence va upsampled con fsampling in PNRcodes 1xN

            PRNsamples=interp1(1:length(obj.PRNsequence),obj.PRNsequence, ...
                1:1/obj.fsampling:length(obj.PRNsequence),"previous"); %upsampling
            
            PRNcode=repmat(PRNsamples,1,windowSize*obj.nPRN_x_Symbol);             
            clear PRNsamples
            
            wsize=size(filteredSamples,1);            
            %parsing
            PRNcode=PRNcode(1:wsize);
            
            %DT$ padding circshift()
            E_PRNcode=[PRNcode(1-round(sampleShifts(1)):end),PRNcode(1:1-round(sampleShifts(1)))];
            %P_PRNcode=PRNcode;
            L_PRNcode=[PRNcode(1:1+round(sampleShifts(3))),PRNcode(1:end-round(sampleShifts(3)))];
            
            %in-phase multicorrelation
            Ie=E_PRNcode.*filteredSamples(:,1)';
            Ip=PRNcode.*filteredSamples(:,1)';
            Il=L_PRNcode.*filteredSamples(:,1)';

            %quadrature multicorrelation
            Qe=E_PRNcode.*filteredSamples(:,2)';
            Qp=PRNcode.*filteredSamples(:,2)';
            Ql=L_PRNcode.*filteredSamples(:,2)';
            
            %parsing
            nCoherentSamples=round(nSamples_x_symbolPeriod*coherenceFraction);
            %roundNSamples=round(size(filteredSamples,1)/nCoherentSamples)*nCoherentSamples;
            
            %DT$
            %corrIQe=[sum(reshape(Ie(1:roundNSamples),nCoherentSamples,[]),1);
            %       sum(reshape(Qe(1:roundNSamples),nCoherentSamples,[]),1)];

            %coherence time reshaping and integration
            Ie=reshape(Ie(1:wsize),nCoherentSamples,[]);
            Ie=sum(Ie,1);
            Ip=reshape(Ip(1:wsize),nCoherentSamples,[]);
            Ip=sum(Ip,1);
            Il=reshape(Il(1:wsize),nCoherentSamples,[]);
            Il=sum(Il,1);

            Qe=reshape(Qe(1:wsize),nCoherentSamples,[]);
            Qe=sum(Qe,1);
            Qp=reshape(Qp(1:wsize),nCoherentSamples,[]);
            Qp=sum(Qp,1);
            Ql=reshape(Ql(1:wsize),nCoherentSamples,[]);
            Ql=sum(Ql,1);

            %non-coherent integration
            %s_IQe=sum(sum(corrIQe.^2,2),1);
            s_IQe=sum(Ie.^2+Qe.^2);
            s_IQp=sum(Ip.^2+Qp.^2);
            s_IQl=sum(Il.^2+Ql.^2);
            %find max
            [~,idShift]=max([s_IQe,s_IQp,s_IQl]); % early -> 1 prompt -> 2 late -> 3            
     
            %decoding
            I=[Ie;Ip;Il];
            %DT$LB
            selI=reshape(I(idShift),1/coherenceFraction,[]); %columns of coherent fractions for each symbol
            decSymbols=(2*(sum(selI,1)>0)-1)'; % column vector of decoded symbols +1,-1            
            
        end

            %somma interna ad un coherence time -> coherent integration
            %somma dei quadrati delle somme risultanti -> non-coherent integration
            %discriminator
            %DOT product: sommatoria[(Ie-Il)*Ip] - sommatoria[(Qe-Ql)*qp],
            %if <0 detector is late, if >0 too early


        function analyzeMessage(obj,decodedSymbols,outData)

            %minimum distance region decision, no channel inversion
            decodedSymbols=int16((decodedSymbols+1)./2); % 1 -> 1   -1 -> 0

            startPoint=1;
            [SyncMessage,startPoint]=extractAndAdvance(decodedSymbols,startPoint,obj.SyncLength);
            [SV_ID,startPoint]=extractAndAdvance(decodedSymbols,startPoint,obj.SV_IDLength);
            [M_ID,startPoint]=extractAndAdvance(decodedSymbols,startPoint,obj.M_IDLength);
            [M_body,startPoint]=extractAndAdvance(decodedSymbols,startPoint,obj.M_bodyLength);
            [CRCMessage,startPoint]=extractAndAdvance(decodedSymbols,startPoint,obj.CRCLength);
            [TailMessage,~]=extractAndAdvance(decodedSymbols,startPoint,obj.TailLength);            
            clear startPoint;

            %M_IDCheck=decodedSymbols(1);
            outData.SV_ID=num2str(SV_ID);
            outData.message_ID=num2str(M_ID);    
            outData.message_Body=num2str(M_body);         
            outData.CRC=num2str(CRCMessage);                      
            
            CRCCheck=checksum([M_ID M_body CRCMessage]);   

            %verifications
            outData.isACKmessage=(M_body(1)==0); 
            outData.ACKed=(sum(CRCCheck)==0);
            %disp("",CRCmessage,CRCCheck)            

            if TailMessage~=obj.Tail 
                disp("error tail"); 
            end

            if SyncMessage~=obj.Sync 
                disp("error sync"); 
            end                       
               
        end


        function CRCCheck = checksum (obj,MessageToCheck)
            %leftmost symbol equal to 1
            %compute checksum only on the superposed portion 
            %CRC
            %G(X) = (X+1)P(X)
            %P(X) = X23 + X17 + X13 + X12 + X11 + X9 + X8 + X7 + X5 + X3 + 1
            %     = 1000 0010 0011 1011 1010 1001
            %     = A23DCB  

            %padding
            tmpMessage=[MessageToCheck,zeros(1,obj.CRCLength)]; %M_ID + M_body + CRC + padding

            j=find(tmpMessage,1); %find symbol 1            
            while j>0 && 1+length(tmpMessage)-j>obj.CRCLength                   
                tmpMessage=tmpMessage(j:end);
                tmpMessage(1:1+obj.CRCLength)=bitxor(tmpMessage(1:1+obj.CRCLength),obj.CRCkey);                
                j=find(tmpMessage,1);                
            end

            CRCCheck=tmpMessage(1+end-obj.CRCLength:end);
            clear tmpMessage

        end


        function [newArray,startPoint]=extractAndAdvance(~,OriginalArray,startPoint,length)
            %prende l'intervallo startPoint-startPoint+length e aggiorna startPoint+=length

            newArray=OriginalArray(startPoint:startPoint+length-1);             
            startPoint=startPoint+length;  


        end
                    
    end
end