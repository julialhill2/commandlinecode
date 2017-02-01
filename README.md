# commandlinecode
scripts for processing EEG/depth electrode data without EEGLAB

To convert data from .edf to .mat format (and for a description of the steps) use "Madineh_edftomat_adapted.m". This is a script that Madineh Sarvestani (a friend from my CSHL class) sent to me. 

"plot_spectrogram.m" is also from Madineh, and will plot spectrograms for several mice. 

"12017_Pwelch_oncontinuoussignal.m" is a script that will allow you to both process the .edf data to .mat, and then do power calculations for a group of files and save the results in an excel file. 

"WORKS_seizures_linelength_FIXED.m" does line length calculations. It bins the selected data, finds the difference between successive data points, sums them and calculates the average value. You have to do this and the sharp waves file for single mice. 

"WORKS_seizures_sharpwaves.m" is a script that I adapted from Henry to identify sharp waves in HPC data. It z-scores the data and then finds points 7 SD above the threshold. You can also see the average length of the events. 
