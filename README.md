# Satellite-lab

# Introduction
Our project centers on the developing an efficient MATLAB simulator for the demodulation of GNSS packets. It is part of the Satellite Communications Systems course and aims at forming on-the-field experience in the simulation of ICT systems and at deepening our understanding of the related problems.

We decide to develop our simulation distancing from the actual FPGA (Field Programmable Gate Array) implementation to exploit the MATLAB capability to parallelize array computations. In particular, we acquire the GNSS signal though a Parallel Code Phase Search (PCPS), and we implement symbol- and segment-wise tracking and decoding. The tracking loop implements a hybrid correlation, i.e., performing incoherent sum of small coherent intervals, to evaluate the peak position at different delays (Code Phase tracking) and Doppler shifts (Doppler Frequency tracking). The demodulation compensates the residual Doppler envelope using a Feed-Forward (FF) Carrier Recovery with a Non-Decision-Aided (NDA) Phase Estimator before performing coherent sum for each symbols of the coherent intervals (Doppler Residual Phase tracking).

Employing together all these techniques, we can evaluate their impact on the complete system and reach performances comparable to the ones of actual GNSS receivers without the complexity of performing real-time (sample-wise) processing.



# IN-OUT INTERFACE
## Input config file

The input setting file is formatted as Json and contains the parameters necessary for the simulation:
  - fSampling: float, 2-60MHz, sampling frequency of the baseband signal (affected by doppler shift)
  - quantizationBits: int, 8/16/24/32, number of bit per I/Q sample in the binary file
  - scenarioDuration: float, 1-60s, duration of the simulated signal (determines the dimension of the IQsamples file)
  - SV_PRN_ID: int, 1-50, PRN sequence of "Primary Codes for the E1-B Component" (see GALILEO OD SIS ICD)
  - CRCpolynomial: hex-string (24bit), generator polynomial G(x) for CRC value creation
  - SYNCpattern: bin-string (10bit), 0101100000, fixed synchronization header
  - TAILpattern: bin-string (6bit) 000000, fixed padding tail
  - SVIDlength: int, 6, length of SV ID field
  - MIDlength: int, 4, length of Message ID field
  - MBODYlength_TX: int, 80, length of Message Body field with TX packets
  - MBODYlength_ACK: int, 30, length of Message Body field with ACK packets
  - CRClength: int, 24, length of CRC sequence
  - nPRN_x_Symbol: int, 1-10, number of PRN repetition for each symbol, determines the symbol period with the following parameters
  - nChip_x_PRN: int, 4092, length of the PRN sequence of "Primary Codes for the E1-B Component" (see GALILEO OD SIS ICD)
  - chipRate: float, 1.023 MHz, chip rate of the PRN sequence
  - maxDoppler: float, 100 KHz, maximum doppler frequency, useful to estimate the signal bandwith and filter noise
  - (maxDopplerShift_x_Symbol: float, 1Hz, maximum doppler shift per symbol period, limit for tracking research)

```
{
  "fSampling": 10e6,
  "quantizationBits": 16,
  "scenarioDuration": 5.05,
  "SV_PRN_ID": 1,
  "CRCpolynomial": "A23DCB",
  "SYNCpattern": "0101100000",
  "TAILpattern": "000000",
  "SVIDlength": 6,
  "MIDlength": 4,
  "MBODYlength_TX": 80,
  "MBODYlength_ACK": 30,
  "CRClength": 24,
  "nPRN_x_Symbol": 1,
  "nChip_x_PRN": 4092,
  "chipRate": 1.023e6,
  "maxDoppler": 1e5
}
```

## Input binary file
The input binary files contains the I/Q samples of the received signals.
Length and number of bits per sample are defined in the config file. Each pair of samples represent the in-phase and quadrature components of the received signal.

```
    00000000      11111111      00000000      11111111      ...
    Isample k=1   Qsample k=1   Isample k=2   Qsample k=2   ...
```

## Output file

The output file is converted directly from a Matlab class into Json format.
It contains all the output parameter from the receiver module:
  - version: string, version of the algorithm
  - ACQUISITION_OK: bool, debug log of acquisition algorithm
  - TRACKING_OK: bool, debug log of the tracking algorithm
  - DEMODULATION_OK: bool, debug log of the message analysis
  - SV_ID: string 6char, binary ID of target satellite (demodulated)
  - message_ID: string 4char, message ID (demodulated)
  - message_body: string 30char, message body (demodulated)
  - CRC: string 24char, CRC of the message (demodulated)
  - ACKed: bool, true if the computed CRC is equal to the received one
  - isACKmessage: bool, true if the message_body contains ACK flag (first bit is 1)
  - estimatedDopplerStart: float, estimated doppler frequency after acquisition procedure
  - estimatedDopplerEnd: float, estimated doppler frequency after message demodulation (tracking procedure)
  - estimatedDelay: float, estimated time delay (from the start of the IQsamples file) during the acquisition procedure
  - estimatedPhase: float, 0-2pi, estimated envelope phase (at estimatedDelay) during the acquisition procedure

```
{
    "version": "2"
    "ACQUISITION_OK", true,
    "TRACKING_OK", true,
    "DEMODULATION_OK", true,
    "SV_ID": "000000",
    "message_ID": "0000",
    "message_body": "000000000000000000000000000000",
    "CRC": "000000000000000000000000",
    "ACKed": true,
    "isACKmessage": false,
    "estimatedDopplerStart": 10.05,
    "estimatedDopplerEnd": 10.05,
    "estimatedDelay": 1.525,
    "estimatedPhase": 0
}
```

# Preliminary Analysis:
## Starting information:
- The input signals are already sampled and acquired (we do not model the acquisition and its effects on the signals) BUT we still need to filter them -> it is equivalent to model pass-band direct acquisition?
- We maybe need to move them into baseband signals to evaluate -> we are sure about what is the input in case the other groups perform erroneous procedures
- The input file is binary, with int8 samples -> we maybe work with integers to save memory space and speed-up computation
- The output is binary sequence (without sync symbols, with CRC)

## Read config file
- It contains the parameters for the receiver
- It must be read from a file and save locally into matlab variables
- The config parameters are input of the processing functions

  ```
  params=struct("version",1,"name","value);
  params=getFromFile("filename"); %you must create the structure CORRECTLY inside the method
  params.version;
  ```

- During the test/development we create local variable and then we asign them from the file when we'll know what it contains

## Read sampling data
- It is a binary file and contains 8bit samples alternated I (in phase) and Q (in quadrature)
- It must be read from a file and save locally IF its dimension is sufficiently small
- The sampling frequency and the time vector is defined from config values

  ```
  IQsamples = zeros(N,2);
  Isamples = IQsamples(:,1);
  Qsamples = IQsamples(:,2);
  t = 0:1/Fs:(N-1)/Fs;
  ```

## Beating with carrier frequency (no doppler)
- f0 = carrier default frequency
- The results has bandwith B (modulation band) + Df (doppler frequency)
- This operation is fictional: should be already performed before acquisition that is not modelled
- We multiply the signal (IQsamples) by cos/sin(2* pi* f0* t)
- TO TEST: verify the equivalence using exp(i* 2* pi* f0* t) (WARNING: superposition of negative spectrum)
- TO TEST: multiply for sqrt(2)
- The output is still a Nx2 signal with the IQ samples

  ```
  IQbeated = beatSignalWithCarrier(IQsamples,carrier);
  ```

## Filter to remove spectrum copies
- FIR/IIR filter (see Matlab passband filters)
- The band of the filter should take into account the expected doppler 

  ```
  IQfiltered = filterSignal(IQbeated);
  ```

## Correlation of the sync portion
- We estimate the sync portion (???)
- We create the signal Ilocal (no doppler, mod data are syncData+PRN) ->+-1 signal with square waveform (or specified in config)
- We calculate the correlation (over a grid of tau x Df) (using FFT and elem-wise multiplication)
- We find the optimal values
- The output are delay and doppler frequency

  ```
  Ilocal = zeros(Nsync,2);
  [tau,Df] =  correlateAndEstimateDelayAndDoppler(IQfiltered,Ilocal);
  ```

## Demodulation
- Using tau and Df we beat the signal with cos/sin(2* pi* Df* t - tau) to obtain baseband signal
- We filter out the spectrum copies (TO TEST)
- We project (inner product) using the PRN (every T=symbol period obtained from Df and tau)
- The output is a sequence Mx2 of projection IQ of symbol points

  ```
  IQrecSymbols = zeros(M,2);
  ```

## Channel inversion
- Using the sync sequence we evaluate the attenuation and phase of the received signal
- We invert the channel

  ```
  IQdemodSymbols = zeros(M,2);
  ```

## Decoding
- We decide with minimum distance the received symbols (+1 and -1)
- We remove the sync sequence head and the CRC
- We compute and verify the CRC
- The output is a binary stream

  ```
  demodBits = zeros(M,1);
  ack = 0;
  ```

## Saving into output file

  ```
  saveBinarySequence("filename",[ack; demodBits]);
  ```


# Implementation choices
- La parte di moltiplicazione e filtraggio la svolgiamo 2 volte perché la prima ha funzione di "validazione" (verificare che sia corretto anche tramite grafici e parametri come la densità spettrale) e pulizia (ridurre rumore out-of-band) dei dati di ingresso, mentre la seconda ci serve per ottenere il segnale modulante e quindi decodificare i simboli.

# Test matrix parameters
- filter bandwidth (B = chip rate):    B,         3B,         no filter
- fragment size (n symbols):           1,         2,          5,          10
- coherent fractions per symbol:       1,         3,          5
- doppler (max + resolution) Hz:       (1e2 + 1), (1e3 + 10), (1e4 + 100)
- delay: 2-3k samples
- doppler start: 45.52 Hz
- amplitudes: noise 2^14 -> signal 2^(14-reduce) [reduce = log2(10.^([power dB]/20))]
- axis X: doppler end-start(symbol) Hz         0(0),     0.2(2.5m),     0.6(7.5m),    1.4(17.5m),      3(37.5m),
                                         6.2(77.5m),  12.6(157.5m),  25.4(317.5m),    51(637.5m), 102.6(1.2825),
                                      132.6(1.6575), 162.6(2.0325),  204.4(2.555),  252.4(3.155),  324.4(4.055),
- axis Y: power signal/noise dB(reduce)   -40(6.64),     -38(6.31),     -36(5.98),     -34(5.65),     -32(5.32),
                                          -30(4.98),     -26(4.32),     -22(3.65),     -18(2.99),     -14(2.32)
                                          -40.5(6.73),   -41(6.81)

We have 108 parameters combinations (3 * 4 * 3 * 3).
We have 100 test scenario (10 * 10).
We generate 15Mb of file for each signal (TOT: 1.5Gb), we need exclude them from GIT but not delete them runtime.

What save at each step: peak value, threshold value, doppler start, doppler end, time delay, demodulation OK, acked message
