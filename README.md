# Passive-Cavitation-Detection
This repository contains the scripts necessary for processing data from passive cavitation detectors.

## Setup
A H107 transducer (0.5 MHz centre frequency) is used to transmit the pressure wave and an unfocused passive cavitation detector (5 MHz centre frequency), coaxially aligned to the transducer,  is used to detect the re-emitted acoustic signal. 

## Files

### periodogram_pcd.m
This file plots a periodogram of the signal. This is to (1) decide what fmin and fmax boundaries should be chosen to calculate the PSD, and (2) decide what frequency range is necessary to capture harmonics and ultraharmonics (if there are any present). 

### call_pcd.m
This file calls the function fprocess_pcd. Here, the parameters below section 1.2 will need to be changed based on settings used and desired output. 

### fprocess_pcd.m
Depending on the results observed by running periodogram_pcd, it might be necessary to change the variables in section 1.4.2 and 1.4.3.


