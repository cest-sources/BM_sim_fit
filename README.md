# BM_sim_fit
Bloch-McConnell simulation (full numerical and analytical) and fit (multi B1)
preRelease 2017

the code was used and is documented for the paper:

Zaiss, M., Angelovski, G., Demetriou, E., McMahon, M. T., Golay, X. and Scheffler, K. (2017), QUESP and QUEST revisited – fast and accurate quantitative CEST experiments. Magn. Reson. Med. doi:10.1002/mrm.26813



Allows for simulation of BM system upon RF presaturation 
 - multiple pools possible (water, 5 CEST, 1  semisolid MT)
 - trains RF pulses with arbitratry pulse shapes possible (gauss, sinc, spin-lock already implemented)
 - both analytical (if possible) and full numerical simulation possible
 - transient and steady-state possible
 - first order readout schemes possible
 - Both full Z-spectra and single offset evaluation possible
 - Ready for multi parameter simulations
 
**Read more about usage in \doc\BM_Documentation.docx 	and \doc\BM_tutorial.pptx**


###Acknowledgments
This code was created by Moritz Zaiss based on the cw BM simulation and fit of Shanrong Zhang. A revised version was created with the help of Patrick Schuenke, Christian Meyer, Christian David and Volkert Roeloffs.
