# Satellite-lab

# Starting information:
- The input signals are already sampled and acquired (we do not model the acquisition and its effects on the signals) BUT we still need to filter them -> it is equivalent to model pass-band direct acquisition?
- We maybe need to move them into baseband signals to evaluate -> we are sure about what is the input in case the other groups perform erroneous procedures
- The input file is binary, with int8 samples -> we maybe work with integers to save memory space and speed-up computation
- The output is binary sequence (without sync symbols, with CRC)


# Processing Steps:
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
- We remove the sync sequence head
- The output is a binary stream

  ```
  demodBits = zeros(M,1);
  ```

## Saving into output file

  ```
  saveBinarySequence("filename",demodBits);
  ```


# Implementation choices
- La parte di moltiplicazione e filtraggio la svolgiamo 2 volte perché la prima ha funzione di "validazione" (verificare che sia corretto anche tramite grafici e parametri come la densità spettrale) e pulizia (ridurre rumore out-of-band) dei dati di ingresso, mentre la seconda ci serve per ottenere il segnale modulante e quindi decodificare i simboli.


