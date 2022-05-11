%
% First implementation: Lorenzo Borsoi
% Review and Testing:
%
classdef Demodulator < handle
    %Demodulator handles...

    properties (SetAccess=public, GetAccess=public)

        CRCkey

        syncLength
        SV_IDLength
        M_IDLength 
        M_bodyLength
        CRCLength 
        TailLength 
        %sync
        %tail

        fsampling
        PNRsequence
        nPRN_x_Symbol

        %windowSize
        %symbolPeriod
        %chipPeriod
        %coherenceFraction
    end

    methods
        function obj = Demodulator()
            %Demodulator constructor
            %   Create DownconverterFilter class
        end

        function obj = configCorrelatorValues(obj,fSampling,nPRN_x_Symbol,PRNsequence)
            obj.fsampling=fSampling;
            obj.nPRN_x_Symbol=nPRN_x_Symbol;
            obj.PNRsequence = PRNsequence;   
        end

        function obj = configMessageAnalyzer(obj,CRCpolynomial, ...
                                            SYNCpattern,SVIDlength,MIDlength, ...
                                            MBODYlength,CRClength,TAILpattern)
            polArray = hexToBinaryVector(CRCpolynomial);
            
            log2dec = @(v) bin2dec(num2str(v));
            sumbin = @(a,b) dec2bin(sum([log2dec(a),log2dec(b)]))-'0';
            
            % P(X) + X*P(X)
            obj.CRCkey = sumbin([0 polArray],[polArray 0]); %test with actual CRC and messages
            
            obj.syncLength=length(SYNCpattern);
            obj.SV_IDLength=SVIDlength;
            obj.M_IDLength=MIDlength; %verify if first bit is zero
            obj.M_bodyLength=MBODYlength;
            obj.CRCLength=CRClength; %compute and verify
            obj.TailLength=length(TAILpattern); %000000
        end

        function [decSymbols,idShift] = decodeOptimumShift(obj,filteredSamples,windowSize, ...
                                                            sampleShifts,nSamples_x_symbolPeriod, ...
                                                            nSamples_x_chipPeriod,coherenceFraction)

            %nchips=length(obj.PNRsequence)*4;
            %PRNcodes=obj.PNRsequence;
            %PRNperiod=length(inout.setting.PRNsequence)/chipPeriod;
            
            %windowSize: numero di simboli
            %wsize: numero campioni
            %nchips: numero di chips che compongono la PNR sequence
            %numero di PNR sequence in un simbolo: symbolPeriod/(chipPeriod*nchips)
            %numero di PNR sequence in un messaggio: windowSize*symbolPeriod/(chipPeriod*nchips)
            %PNRsequence va upsampled con fsampling in PNRcodes 1xN

            %k=fsampling*chipPeriod=10.2
            %mod(k,1)=0.2 %numero dopo la virgola
            %approssimazione alla quarta cifra?
            %1/0.2=5 per ogni ciclo di chip, il quinto ha un sample in più
            %tronco al numero intero
            %repmat(PNRcodes,k_int,1);
            %aggiungo ciclicamente il sample in più al quinto chip
            %reshape
            % 10 10 10 10 11
            
            %DT resample oppure upsample/downsample  100 e int(chipPeriod/100)
            %   repelem(PRNcodes,round(chipPeriod))
            
            PRNcodes=repmat(obj.PNRsequence,1,windowSize*obj.nPRN_x_Symbol); 
            
            %upsampling (non funziona)
            %[p,q] = rat(obj.fsampling / (1/obj.chipPeriod)); %manca chipPeriod, da controllare
            [p,q] = rat(nSamples_x_chipPeriod);
            PRNcodes=downsample(upsample(PRNcodes,p),q); 

            %TODO
            wsize=size(filteredSamples,1);            
            %parsing
            PRNcodes=PRNcodes(1:wsize);
            
            %padding circshift()
            %E_PNRcodes=[PNRcodes(round(sampleShifts(3)+1):end),zeros(round(sampleShifts(3)),1)];
            E_PNRcodes=[PNRcodes(1-round(sampleShifts(1)):end),PNRcodes(1:1-round(sampleShifts(1)))];
            %P_PRNcodes=PRNcodes;
            L_PNRcodes=[PNRcodes(1:1+round(sampleShifts(3))),PNRcodes(1:end-round(sampleShifts(3)))];
            
            %in-phase multicorrelation
            Ie=E_PNRcodes.*filteredSamples(:,1)';
            Ip=PRNcodes.*filteredSamples(:,1)';
            Il=L_PNRcodes.*filteredSamples(:,1)';

            %quadrature multicorrelation
            Qe=E_PNRcodes.*filteredSamples(:,2)';
            Qp=PRNcodes.*filteredSamples(:,2)';
            Ql=L_PNRcodes.*filteredSamples(:,2)';
            
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

            %da riscrivere così
                %decodedSymbols(1:obj.syncLength)
                %decodedSymbols(1+obj.syncLength:obj.syncLength+obj.SV_IDLength)
                %startPoint=1;
                %[x,startPoint]=extractAndAdvance(startPoint,obj.syncLength);%prende l'intervallo startPoint-startPoint+length e aggiorna startPoint+=length
                %y=extractAndAdvance(startPoint,obj.SV_IDLength);

            TailMessage=decodedSymbols(1+end-obj.TailLength:end);
            SyncMessage=decodedSymbols(1:obj.syncLength);
            decodedSymbols=decodedSymbols(1+obj.syncLength:end-obj.TailLength); %remove sync remove tail
            
            outData.SV_ID=num2str(decodedSymbols(1:obj.SV_IDLength));
            decodedSymbols=decodedSymbols(1+obj.SV_IDLength:end); %remove SV_ID
            
            CRCMessage=decodedSymbols(1+end-obj.CRCLength:end);             
            decodedSymbols=decodedSymbols(1:end-obj.CRCLength); %M_ID + M_body 
            %M_IDCheck=decodedSymbols(1);

            outData.message_ID=num2str(decodedSymbols(1:obj.M_IDLength));          
            outData.message_Body=num2str(decodedSymbols(1+obj.M_IDLength:end));           
            outData.CRC=num2str(CRCMessage);
   
            %CRC
            %G(X) = (X+1)P(X)
            %P(X) = X23 + X17 + X13 + X12 + X11 + X9 + X8 + X7 + X5 + X3 + 1
            %     = 1000 0010 0011 1011 1010 1001
            %     = A23DCB            
            
            CRCCheck=checksum(decodedSymbols);   

            %verifications
            if sum(CRCCheck)==0
                outData.ACKed=true;      
            else
                outData.ACKed=false; 
                %disp("",CRCmessage,CRCCheck)
            end

            if TailMessage~=obj.tail 
                disp("error tail"); 
            end %(da creare in config)
            if SyncMessage~=obj.sync 
                disp("error sync"); 
            end 

            outData.isACKmessage=(decodedSymbols(1+obj.M_IDLength)==0); % $LB$ sicuro?            
               
        end

        function CRCCheck = checksum (decodedSymbols)

            %leftmost symbol equal to 1
            %compute checksum only on the superposed portion 

            %padding
            tmpMessage=[decodedSymbols,zeros(1,obj.CRCLength)]; %M_ID + M_body + CRC + padding

            j=find(tmpMessage,1); %find symbol 1
            
            while j>0 && 1+length(tmpMessage)-j>obj.CRCLength                   
                tmpMessage=tmpMessage(j:end);
                tmpMessage(1:1+obj.CRCLength)=bitxor(tmpMessage(1:1+obj.CRCLength),obj.CRCkey);                
                j=find(tmpMessage,1);                
            end

            CRCCheck=tmpMessage(1+end-obj.CRCLength:end);

        end
                    
    end
end