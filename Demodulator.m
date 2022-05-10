%
% First implementation: Lorenzo Borsoi
% Review and Testing:
%
classdef Demodulator < handle
    %Demodulator handles...

    properties (SetAccess=public, GetAccess=public)
        symbolPeriod
        PNRsequence
        fsampling
        chipPeriod
        coherenceFraction
    end

    methods
        function obj = Demodulator()
            %Demodulator constructor
            %   Create DownconverterFilter class
        end

        function obj = configCorrelatorValues()
            
        end

        function obj = configMessageAnalyzer(key)
            obj.PNRsequence=inout.setting.PRNsequence;    
            obj.key=key;
        end

        function [decSymbols,idShift] = decodeOptimumShift(filteredSamples,sampleShifts)

            %nchips=length(obj.PNRsequence)*4;
            PRNcodes=obj.PNRsequence;
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

            PRNcodes=repmat(PRNcodes,1,windowSize*obj.symbolPeriod/(obj.chipPeriod*length(obj.PNRsequence)*4));
            %......
            wsize=size(filteredSamples(:,1),1);
            %parsing
            PRNcodes=PRNcodes(1:wsize);
            
            %padding
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
            coherentSamples=round(wsize/windowSize/obj.coherentFraction);
            %coherence time reshaping and integration
            Ie=Reshape(Ie(1:wsize),wsize/windowSize/obj.coherentFraction,[]);
            Ie=sum(Ie,1);
            Ip=Reshape(Ip(1:wsize),wsize/windowSize/obj.coherentFraction,[]);
            Ip=sum(Ip,1);
            Il=Reshape(Il(1:wsize),wsize/windowSize/obj.coherentFraction,[]);
            Il=sum(Il,1);

            Qe=Reshape(Qe(1:wsize),wsize/windowSize/obj.coherentFraction,[]);
            Qe=sum(Qe,1);
            Qp=Reshape(Qp(1:wsize),wsize/windowSize/obj.coherentFraction,[]);
            Qp=sum(Qp,1);
            Ql=Reshape(Ql(1:wsize),wsize/windowSize/obj.coherentFraction,[]);
            Ql=sum(Ql,1);

            %non-coherent integration
            s_IQe=sum(Ie.^2+Qe.^2);
            s_IQp=sum(Ip.^2+Qp.^2);
            s_IQl=sum(Il.^2+Ql.^2);
            %find max
            [~,idShift]=max([s_IQe,s_IQp,s_IQl]); % early -> 1 prompt -> 2 late -> 3
            
            %take the sum of the corrValues in the currentSymbol time period
            %decoding
            I=[Ie;Ip;Il];
            if sum(I(idShift,(currentSymbol-1)*coherentFraction+1:currentSymbol*coherentFraction))>=0
                decSymbol=1;
            else 
                decSymbol=-1;
            end            
            
        end

            %somma interna ad un coherence time -> coherent integration
            %somma dei quadrati delle somme risultanti -> non-coherent integration
            %discriminator 
            %DOT product: sommatoria[(Ie-Il)*Ip] - sommatoria[(Qe-Ql)*qp],
            %if <0 detector is late, if >0 too early


        function analyzeMessage(decodedSymbols)

            %minimum distance region decision, no channel inversion
            decodedSymbols=(decodedSymbols+1)./2;
            %get struct
            syncLength=10;
            SV_IDLength=6;
            M_IDLength=4; %verify if first bit is zero
            M_bodyLength=30;
            CRCLength=24; %compute and verify
            TailLength=6; %000000

            TailMessage=decodedSymbols(1+end-TailLength:end);
            decodedSymbols=decodedSymbols(1+syncLength:end-TailLength); %remove sync remove tail
            
            outData.SV_ID=num2str(decodedSymbols(1:SV_IDLength));
            decodedSymbols=decodedSymbols(1+SV_IDLength:end); %remove SV_ID
            
            CRCMessage=decodedSymbols(1+end-CRCLength:end);             
            decodedSymbols=decodedSymbols(1:end-CRCLength); %M_ID + M_body 
            M_IDCheck=decodedSymbols(1);
                        
            %CRC
            %inout.setting.CRCCheck

            %leftmost symbol equal to 1
            %compute checksum only on the superposed portion 
            
            %G8X) = (X+1)P(X)
            %P(X) = X23 + X17 + X13 + X12 + X11 + X9 + X8 + X7 + X5 + X3 + 1
            %padding
            newMessage=[decodedSymbols,zeros(1,CRCLength)]; %M_ID + M_body + CRC + padding

            j=find(newMessage,1); %find symbol 1
            
            while j>0 && 1+length(newMessage)-j>CRCLength                   
                newMessage=newMessage(j:end);
                newMessage(1:1+CRCLength)=bitxor(newMessage(1:1+CRCLength),obj.key);                
                j=find(newMessage,1);                
            end

            CRCCheck=newMessage(1+end-CRCLength:end);
           
            %verification
            if (M_IDCheck==0 && CRCCheck==zeros(1,CRCLength) && TailMessage==zeros(1,TailLength))
                outData.ACKed=true;            
            else
                outData.ACKed=false;          
            end         
                       
            
            outData.message_ID=num2str(decodedSymbols(1:M_IDLength));          
            outData.message_Body=num2str(decodedSymbols(1+M_IDLength:end));           
            outData.CRC=num2str(CRCmessage);
            %outData.CRC=num2str(CRCCheck);
               
        end
    end
end