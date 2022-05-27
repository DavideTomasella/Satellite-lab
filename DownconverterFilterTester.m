%% PARAMETERS
clearvars
close all

inout = InOutInterface();

% directories
indir = "inData";
inout.configCreateSettings(indir);
outdir = "outData";
inout.configSaveResults(outdir);

filename = "in0.json";
PRNfilename = "PRNpattern.json";
settings = inout.createSettings(filename,PRNfilename);

down = DownconverterFilter();
down.configDownConverter(inout.settings.fSampling);

signal = SignalManager();
signal.configReadFile("binData","nine.bin",inout.settings.quantizationBits);

ampl = 2^(inout.settings.quantizationBits-2);   % gain of the signal, otherwise 
                                                % the in16 matrix results made of 1,0 and -1
dop = 10E3;         % doppler frequency 10E3
sigma = 1E4;        % noise variance 1000
initial_delay = 0;  %initial delay in seconds

%% PRN SEQUENCE GENERATION
tchip = (0:length(inout.PRNcode)-1)' / settings.chipRate;
upsampling = settings.fSampling / settings.chipRate;
t = (0:upsampling*length(inout.PRNcode)-1)' / settings.fSampling;
PRN = -ampl*interp1(tchip,2*inout.PRNcode-1,t,"previous");
PRN = PRN(1:end-round(upsampling));
t = t(1:end-round(upsampling));

%% BASEBAND SIGNAL GENERATION WITH NOISE AND DOPPLER
% delay = sigma*randn(1,initial_delay*settings.fSampling)';

IQ = testSignal(t,PRN,dop,initial_delay,sigma);

figure;
plot(t,IQ(:,1));
hold on;
plot(t,PRN);
title("PRN and modulated signal (baseband) with added noise");
xlim([0.0013 0.0014]);

signal.IQsamples_float = IQ;

%% USING A TEST SIGNAL
signal.configReadFile("binData\testSignals","T_tracking_24b.bin",inout.settings.quantizationBits);
signal.readFile(1,400000);
t = (0:length(signal.IQsamples)-1)/inout.settings.fSampling;
dop = 418.7;

%% FILTERING THE BASEBAND SIGNAL
down.configFilter(1,250);
down.downFilter(signal,0.6*inout.settings.chipRate,inout.settings.chipRate);

figure;
subplot(2,1,1);
plot(t,signal.IQsamples(:,1));
hold on;
% plot(t,k(1:length(reader.IQsamples)));
title("PRN and filtered baseband signal");
subplot(2,1,2);
plot(t,signal.IQsamples(:,1));
xlim([0.0013 0.0014]);

%% DOWN CONVERSION AND FILTERING
phase = 0;
down.downConverter(signal,dop,initial_delay,phase);

down.configFilter(1,250);
down.downFilter(signal,0.6*inout.settings.chipRate,inout.settings.chipRate);

figure;
subplot(2,1,1);
plot(t,signal.IQsamples_float(:,1)+ampl);
hold on;
% plot(t,k-ampl);
title("Initial PRN and down converted sequence");
subplot(2,1,2);
plot(t,signal.IQsamples_float(:,1)+ampl);
xlim([0.0013 0.0014]);

%% TEST WITH DIFFERENT ERRORS IN DELAY AND DOPPLER
phase = 0;
inital_doppler = 10E3;
doppler_step = 100;
delay_step = 1/(10*inital_doppler);
doppler_sweep = (0:5)*doppler_step + inital_doppler;
delay_sweep = (0:5)*delay_step + initial_delay;
signal.IQsamples_float = IQ;
t = (0:upsampling*length(inout.PRNcode)-1)' / settings.fSampling;
t = t(1:end-round(upsampling));
% sweep in frequency
figure;
plot(t,PRN,"DisplayName","PRN","LineWidth",2);
hold on;
for f_doppler = doppler_sweep
    test_downConverterFilter(t,signal,f_doppler,delay_sweep(1), ...
                                inout.settings.fSampling,inout.settings.chipRate, ...
                                        0,inital_doppler,initial_delay,phase);
    signal.IQsamples_float = IQ;
end
xlim([0.0013 0.0014]);
legend;

% sweep in delay
figure;
plot(t,PRN,"DisplayName","PRN","LineWidth",2);
hold on;
for t_delay = delay_sweep
    test_downConverterFilter(t,signal,doppler_sweep(1),t_delay, ...
                                inout.settings.fSampling,inout.settings.chipRate, ...
                                        1,inital_doppler,initial_delay,phase);
    signal.IQsamples_float = IQ;
end
legend;
xlim([0.0013 0.0014]);

%% FILTER TEST
% B=0.6*chipFrequency & att=250dB/dec is good
% B=0.7*chipFrequency & att=400dB/dec reduce too much harmonics
% B=0.6*chipFrequency & att=350dB/dec reduce too much harmonics (near
% fundamental)
passband = 1*inout.settings.chipRate;
input = testSignal(t,PRN,0,0,sigma);

testFilter = DownconverterFilter();
testFilter = testFilter.configFilter(1,100,inout.settings.fSampling);

signal.IQsamples_float = input;
signal = testFilter.downFilter(signal,passband,inout.settings.chipRate);
Out_with_noise = signal.IQsamples_float;

signal.IQsamples_float = [PRN PRN];
signal = testFilter.downFilter(signal,passband,inout.settings.chipRate);
Out = signal.IQsamples_float;

figure;
subplot(3,1,1);
plot(t,PRN);
hold on;
plot(t,Out(:,1));
xlim([0.001305 0.00132]);
legend("initial","initial filtered",Location="south");

subplot(3,1,2);
plot(t,input(:,1));
hold on;
plot(t,Out_with_noise(:,1));
xlim([0.001305 0.00132]);
legend("initial with noise","initial with noise filtered",Location="south");
subplot(3,1,3);
plot(t,input(:,1));
hold on;
plot(t,Out_with_noise(:,1));
xlim([0.001305 0.00132]);
ylim([0 4E4])
legend("initial with noise","initial with noise filtered",Location="south");

%% SNR PERFORMANCES
MaxDop = 50E3;
fresolution = 1000;
sigma = 1E4;
t = (0:length(PRN)-1)'/inout.settings.fSampling;
f = (0:fresolution:MaxDop);
att = (50:50:450);
testFilter = DownconverterFilter();
% testFilter = testFilter.configFilter(1,350,inout.settings.fSampling);
SNR_0 = zeros(1,length(f));
SNR_F = zeros(length(att),length(f));
SSR = zeros(length(att),length(f));

for h = 1:length(f)
    passBand = 1*inout.settings.chipRate+f(h);
    sig = testSignal(t,PRN,f(h),0,0);
    noise = sigma*randn(length(t),2);
    SNR_0(h) = 10*log10(sum(sig(:,1).^2))-10*log10(sum(noise(:,1).^2));
    
    for g = 1:length(att)
        disp("PassBand = "+num2str(passBand)+" Attenuation = " + num2str(att(g)))
        testFilter = testFilter.configFilter(1,att(g),inout.settings.fSampling);
        signal.IQsamples_float = sig;
        signal = testFilter.downFilter(signal,passBand,inout.settings.chipRate);
        sigFilt = signal.IQsamples_float;
        signal.IQsamples_float = noise;
        signal = testFilter.downFilter(signal,passBand,inout.settings.chipRate);
        noiseFilt = signal.IQsamples_float;
        SNR_F(g,h) = 10*log10(sum(sigFilt(:,1).^2))-10*log10(sum(noiseFilt(:,1).^2));
        SSR(g,h) = 10*log10(sum(sigFilt(:,1).^2))-10*log10(sum(sig(:,1).^2));
%         if mod(h,10) == 0
%             figure;
%             plot(t,sig(:,1)+noise(:,1));
%             hold on;
%             plot(t,sigFilt(:,1)+noiseFilt(:,1));
% %             plot(t,signal.IQsamples_float(:,1)+1E4);
%             grid on;
%             xlim([0.0013 0.0014]);
%             title("Doppler Frequency "+num2str(f(h))+ "Attenuation " + num2str(att(g)));
%         end
    end
end

figure;
plot(f/1E3,SNR_0,"DisplayName","Without Filter");
hold on;
for g = 1:length(att)
    figure;
    plot(f/1E3,SNR_F(g,:),"DisplayName",num2str(-att(g))+" dB/dec");
    legend;
end
legend;
xlabel("Doppler Frequency [kHz]");
ylabel("SNR [dB]");
title("Signal to Noise Ratio");
grid on;

figure;
yline(0,'--r',"DisplayName","Without Filter");
hold on;
for g = 1:length(att)
    plot(f/1E3,SSR(g,:),"DisplayName",num2str(-att(g))+" dB/dec");
end
legend;
xlabel("Doppler Frequency [kHz]");
ylabel("Attenuation [dB]");
title("Signal attenuation");
ylim([-2 0.5]);
grid on;

%% PHASE TEST
fdop = 10E3;
IQ = testSignal(t,PRN,fdop,0,0);
phase = (0:pi/3:2*pi);

figure;
plot(t,PRN,"DisplayName","Initial");
hold on;
for k = (1:length(phase))
    signal.IQsamples_float = IQ;
    signal = down.downConverter(signal,fdop,0,phase(k));
    % the down converted signal has a sign given by the cos(phi) for the
    % In-Phase component and by sin(phi) for the Quadrature component
    plot(t,circshift(sign(cos(phase(k)))*signal.IQsamples_float(:,1),k),"DisplayName",num2str(phase(k)/pi)+" \pi");
end
grid on;
xlim([0.0013 0.0014]);
legend;

%% LOCAL FUNCTIONS
function IQ = testSignal(t,PRN,doppler,delay,sigma)
    noise = sigma*randn(1,length(PRN))';
    basebandsignaldopp = (PRN+noise).*exp(1i*2*pi*doppler*(t-delay));
    IQ = [real(basebandsignaldopp) imag(basebandsignaldopp)];
end