# Satellite-lab

## Starting information:
- The input signals are already sampled and acquired (we do not model the acquisition and its effects on the signals) BUT we still need to filter them -> it is equivalent to model pass-band direct acquisition?
- We maybe need to move them into baseband signals to evaluate -> we are sure about what is the input in case the other groups perform erroneous procedures
- The input file is binary, with int8 samples -> we maybe work with integers to save memory space and speed-up computation
- The output is binary sequence (without sync symbols, with CRC)


## Processing Steps:
### Read sampling data
