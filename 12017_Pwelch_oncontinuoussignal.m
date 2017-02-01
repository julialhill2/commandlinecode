% Converting .edf file to .mat for power spectra and spectrogram
clear all
clc
mainpath='D:\JH\ephys\10.16_Combined_allCSTTrkBfiles_bydate\P30\';
filepath='D:\JH\ephys\10.16_Combined_allCSTTrkBfiles_bydate\P30\';
cd(mainpath);
addpath(genpath(mainpath));


files=dir([filepath, '*.edf']); %find all the .edf files in folder
nfiles=length(files); %number of files in folder
cutOffFreq = 55; % Hz
filterOrder = 3; % Filter order (e.g., 2 for a second-order Butterworth filter)
FS=2000; %even though this is the same for all files/channels, it's good practice to carry it with data
[b, a] = butter(filterOrder, cutOffFreq/(FS/2)); % Generate filter coefficients

%[b,a]=butter(n,Wn)
%returns the transfer function coefficients of an nth-order lowpass digital
%Butterworth filter with normalized cutoff frequency Wn

for i=1:nfiles
    
    mouse_id{i}=files(i).name(1:6); %pull out the 4 letter identifier- can change based on the name length

    [hdr, record] = edfread(files(i).name); %pull out data
    
    HPC = filtfilt(b, a, record(1,:)); % Apply filter to raw EEG1
    EEG2 = filtfilt(b, a, record(2,:)); % Apply filter to raw EEG1
    EMG = filtfilt(b, a, record(3,:)); % Apply filter to raw EEG2

    clear record
    
    filepath='D:\JH\ephys\9.16 ERP_Tcf4\RAWDATA_ALLGROUPS\SCN_allgroups\'; % a subfolder inside the main fodler that holds the .edf (and produced) .mat files
    mousename=mouse_id{i};
    %save the contents to a .mat file- variables listed in quotes will be
    %retained for each animal
    save([filepath,files(i).name(1:6),'.mat'],'EEG1','EEG2', 'EMG', 'FS', 'mousename')
    clearvars -except files nfiles cutoffFreq filterOrder b a FS 
end

%%
clear all
clc
mainpath='D:\JH\ephys\10.16_Combined_allCSTTrkBfiles_bydate\P32\';
filepath='D:\JH\ephys\10.16_Combined_allCSTTrkBfiles_bydate\P32\';
cd(mainpath);
addpath(genpath(mainpath));


files=dir([filepath,'*.mat']); %find all the .mat files in folder
nfiles=length(files);
FS=2000;

for i=1:nfiles
load(files(i).name);
    
start_time=0.1;
stop_time=210;
%In minutes, change time to only include data from a specific time period-
%has to be consistent across channels and across individuals

HPC=single(HPC);
EEG1=EEG1((start_time*60*FS):(stop_time*60*FS));


window=FS*1;
noverlap=floor(0.5*window);
nfft=2^nextpow2(window);

[pxx,f] = pwelch(HPC,window,noverlap,nfft,FS);

spectra = (10*log10(pxx));

%plot(f,spectra); xlim([0 50])

%This exports the dB scaled results to an excel file. The file name will
%save at the top, and the frequencies as well as the spectra will be in
%there. Each sheet is 1 animal. 
filename= '1.31.17_P32HPCpowerspectrapart1';
filename = [filename, '.xlsx'];
nsubj= '14';
xlswrite(filename, mousename, i, 'A1')
xlswrite(filename, f, i, 'A2')
xlswrite(filename, spectra, i, 'B2')
xlswrite(filename, length(HPC), i, 'C2')


        
end
