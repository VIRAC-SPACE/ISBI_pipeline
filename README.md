# Table of Contents
- [Capabilities of ISBI_pipeline](#capabilities)
- [Dependencies](#dependencies)
- [Changelog](#changelog)
- [Getting Help](#getting-help)
- [Acknowledgements](#acknowledgements)
  
## Capabilities
The line and continuum data sets are first loaded. Beginning with the continuum data, corrections for the losses during digital sampling are applied. Then a-priori gain calibration tables (system temperature and gain curve) derived from noise diode temperature measurements conducted during observations are used to calibrate the flux density of visibilities. Bandpass corrections for all targets are made based on the observed bandpass shapes of continuum calibrators. Then, three stages of fringe fitting and integrations are performed. Starting with a `manual phase-cal' stage, the phase difference between the RCP and LCP data are corrected, and any phase delay difference between the baseband channels is corrected allowing the channels and both polarisations to be integrated to improve signal to noise. The group delay is then determined on continuum sources by fringe fitting on the partially integrated data at a solution interval of about 20 minutes in order to trace the slowly drifting delay differences and slowly changing phase residuals imparted by the ionosphere. Solutions are applied to all targets in the experiment. Finally the continuum sources are then fringe fitted again and in so dealing with the baseband channels individually to obtain channels specific solutions. At this stage the absolute flux scale for the experiment is fine tuned by comparison of the measured fluxes of maser and continuum calibrators, providing corrections where needed. It should be noted that all fringe fitting stages up to this point are instructed to discard phase rate solutions. 

At this point, the solutions for the continuum data baseband channel which matches in frequency to the line data set is copied to the line data, thus providing gain, delay and slowly changing phase solutions which can be interpolated to the timeranges of the maser target data. The peak channel of maser emission for each source is then determined and used as the input of a fringe fitting stage which determines the phase and rate fluctuations of the atmosphere with 20 second solution intervals. Phase and rate solutions are concatenated into a single solution table and copied back to the continuum data set, thus enabling long integrations aiming to detect the (typically weak) continuum emission associated with maser targets. A final fringe fitting stage is conducted on the continuum data of all targets as an inspection step as the success or failure at this stage indicates the overall detectability of continuum emission in all sources, both quasars and high-mass protostars. Finally the spectra and visibility plots of all sources are output and the integrated radio continuum flux densities of all sources determined.

## Dependencies
- AIPS
- ParselTongue


## Changelog
Version 1.0

# Getting Help

Bug reports, feature requests and make contributions (e.g. code patches) can be reported by opening a &quot;new issue&quot; ticket on GitHub. Please give as much information (e.g. the software component, version) as you can in the ticket. For bugs, it is extremely useful if a small self-contained code snippet that reproduces the problem is provided.

## Acknowledgements
This software was written by Ross Alexander Burns. If you make use of this software to get results that appear in a publication or presentation please include this acknowledgement: &quot;We have made use of ISBI_pipeline
, a tool developed by Ross Alexander Burns.&quot;

This work was supported by Latvian Council of Science Project "A single-baseline radio interferometer in a new age of transient astrophysics" Nr.: lzp-2022/1-0083.
