# Bleach Encoding

[![DOI](https://zenodo.org/badge/928951568.svg)](https://doi.org/10.5281/zenodo.17520158)

Matlab class for simulating and evaluating how to encode data into fluorescently labeled filaments by photobleaching.

## Installation

Clone this repository:

```bash
git clone https://github.com/thawn/bleach-encoding.git
```

In order to be able to open microscopy images, you need to download and extract the [bioformats MATLAB Toolbox](https://www.openmicroscopy.org/bio-formats/downloads/) into the base directory of this repository.


## Getting started:

the following MATLAB code will generate a simulated microtubule, bleach a pattern and analyze the result:

```matlab
BES = BleachEncodingSim();
%plot a linescan before bleaching
BES.plotLineScan();
%bleach 11 maxima into the microtubule for 0.2s
BES.bleachWaveCounts(11, 0.2);
% plot the linescan after bleaching
BES.plotLineScan();
% plot the result of the fast fourier transform analysis
BES.plotFFTAnalysisWaveNumbers();
```
