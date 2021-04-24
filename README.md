# Mutli_CVS_model
Multiscale Model of the Human Cardiovascular System

This repository contains the MATLAB code to run simulations with a multiscale model of the human cardiovascular system. More informations on the model avalaible here: https://doi.org/10.1016/j.mbs.2016.05.007

The repository is organized as follows:

* Baseline_CVS - contains the code to run the baseline code.
* IIP
  * IIP - contains the code to run Instantaneous Increase in Preload (IIP) protocols, with and without Length-Dependent Activation (LDA)
  * FS_curve - contains the code to run a Frank-Starling curve
* Vascular filling - contains the code to run vascular filling simulations, with and without LDA

Any simulation starts by running the "RunAndSensitivity.m" MATLAB file. 
