# SNR_as_function_B0
This set of scripts have been created to generate the figures in the "Low Field MRI: an MR Physics perspective" published in JMRI in 2019.

In there it is shown that the power law best describing the SNR field dependence in MRI is contrast type dependent (T1, PD and T2*) as well as sequence acquisition strategy (2D vs 3D) dependent. The effects of increasing longitudinal relaxation times (T1) and decreasing apparent transverse relaxation times (T2*) effectively decreasing the SNR gain expected.

Run the script SNRandCNRfieldDependceJMRI.m from its directory.

Currently the SNR and CNR dependence as a function of field strength are estimated assuming the relaxation times of grey and white matter as modelled by (Rooney et al MRM 2007) and (Pohmann et al MRM 2016).
The assumed Power Law of the NMR SNR is 3/2. This assumes the detectable signal increases with the square power of the static magnetic field and noise is in a regime in between system and physiological dominated ones.


The sequence acquisition that is assumed in this simulations assumes

|<----    dead time, dt    ---->|<----    acquisition time,  1/BW    ---->|

|<----    repetition time, TR, divided by number of slices           ---->|

|<----  dt/2 ---->|<----    echo time , TE      -->|<----   1/BW/2   ---->|

"dead time" is used here in the loose sense, means simply that the that acquisition is not taking place (crusher gradients, slice selective gradient, phase encoding, readout preaphaser gradients arebeing  aplpied in this time). The centre of the excitation pulse happens in the middle of this "dead time".  
