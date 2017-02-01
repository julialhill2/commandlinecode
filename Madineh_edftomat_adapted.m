%Script to convert .edf files to .mat data, filter .mat data, and save it to the same directory
%containing the .edf files.

%This script should be called before edftomat.m, which uses the output .mat
%files of this script.

%First part uses the function (edfread.m) to convert .edf files to .mat files for a sample file.
    %Then the sample data is filtered, and nonrelevant channels are deleted.
%The second part creates a big for loop to convert all .edf files to .mat
%files, clean them up, and save them in the same directory.

%Usage:
%To use this script, you have to to have a main directory (called mainpath
%below) that holds all of the .m files. Inside this directory, there should
%be a folder for the raw .edf files, which will be used to store the
%converted .mat files.

%Change lines 29-30 below to set the mainpath and the filepath for your PC.

%Madineh Sarvestani 2015 : continuation of discussion at CSHL Neural Data
%Science

%% Setup Matlab
clearvars
close all
clc
warning off %turns off annoying warning messages
%% Set paths
mainpath='D:\JH\ephys\10.16_Combined_allCSTTrkBfiles_bydate\P32\'; %update this with your path
%mainpath='C:\Users\Madineh\Dropbox\Julie\'; %update this with your path
filepath='D:\JH\ephys\10.16_Combined_allCSTTrkBfiles_bydate\P32\'; % a subfolder inside the main fodler that holds the .edf (and produced) .mat files
%must have a backslash for the folder otherwise you will get error 'index
%exceeds matrix dimensions'

cd(mainpath);filepath
addpath(genpath(mainpath));
%% Let's start wth a sample file

files=dir([filepath, '*.edf']); %find all the .edf files in folder
nfiles=length(files); %number of files in folder

for i=1:nfiles
    %:nfiles
%for i=1; %for just 1 file
[hdr, record] = edfread(files(i).name); %pull out the records
end

%if this doesn't work, but there is a cell for the file, can put
%files(cell#).name)

FS=hdr.samples(1); %pull out sampling rate

%% Setting up the FFT variables
EEG1=record(1,:);
EEG2=record(2,:);
HPC=record(3,:);
%record(3,:);

L=length(EEG1);
NFFT = 2^nextpow2(L); % number of points for NFFT: Next power of 2 from length of y
EEG1_fft = fft(EEG1,NFFT)/L; %fft the EEG1 signal
EEG2_fft = fft(EEG2,NFFT)/L;
HPC_fft = fft(HPC,NFFT)/L;

f = FS/2*linspace(0,1,NFFT/2+1); %setup frequency

%% Filter the data 2: Setup filters and do filtering

cutOffFreq = 50; % Hz
filterOrder = 3; % Filter order (e.g., 2 for a second-order Butterworth filter)
[b, a] = butter(filterOrder, cutOffFreq/(FS/2)); % Generate filter coefficients
filteredEEG1 = filtfilt(b, a, EEG); % Apply filter to EEG1
filteredEEG2 = filtfilt(b, a, HPC); % Apply filter to EEG2

EEGfilt_fft = fft(filteredEEG1,NFFT)/L; %fft the EEG1 signal
HPCfilt_fft = fft(filteredEEG2,NFFT)/L;
%%
figure;
a(1)=subplot(3,2,1);
plot(f,2*abs(EEG_fft(1:NFFT/2+1))); title('EEG1 sepctral');
a(2)=subplot(3,2,3);
plot(f,2*abs(HPC_fft(1:NFFT/2+1))); title('HPC spectra');
%a(3)=subplot(3,2,5);
%plot(f,2*abs(EMG_fft(1:NFFT/2+1))); title('EMG spetra');
xlabel('Frequency (Hz)');
linkaxes([a],'x'); %lock the x-axis together for subplots so you can zoom all
xlim([0 100]); %only look up to 100 Hz

b(1)=subplot(3,2,2);
plot(f,2*abs(EEGfilt_fft(1:NFFT/2+1))); title('Filtered EEG1 sepctral');
b(2)=subplot(3,2,4);
plot(f,2*abs(HPCfilt_fft(1:NFFT/2+1))); title('Filtered HPC spectra');
%b(3)=subplot(3,2,6);
%plot(f,2*abs(EMGfilt_fft(1:NFFT/2+1))); title('Filtered EMG spetra');
xlabel('Frequency (Hz)');
linkaxes([b],'x'); %lock the x-axis together for subplots so you can zoom all
xlim([0 100]); %only look up to 100 Hz

%% Now pass back the cleaned up data for further processing
% Now that we know what the data looks like and how to filter, we'll change
% all the .edf files to cleaned up .mat files in one big for loop without
% any figure output
%clearvars -except filepath mainpath %get rid of everything except the path variables
%clc
%close all

mainpath='D:\JH\ephys\10.16_Combined_allCSTTrkBfiles_bydate\P32\'; %update this with your path
%mainpath='C:\Users\Madineh\Dropbox\Julie\'; %update this with your path
filepath='D:\JH\ephys\10.16_Combined_allCSTTrkBfiles_bydate\P32\'; % a subfolder inside the main fodler that holds the .edf (and produced) .mat files
%must have a backslash for the folder otherwise you will get error 'index
%exceeds matrix dimensions'

FS=2000
files=dir([filepath, '*.edf']); %find all the .edf files in folder
nfiles=length(files); %number of files in folder
%cutOffFreq = 100; % Hz
%filterOrder = 3; % Filter order (e.g., 2 for a second-order Butterworth filter)
%[b, a] = butter(filterOrder, cutOffFreq/(FS/2)); % Generate filter coefficients

%[b,a]=butter(n,Wn)
%returns the transfer function coefficients of an nth-order lowpass digital
%Butterworth filter with normalized cutoff frequency Wn

for i=13
    %:nfiles
    
    mouse_id{i}=files(i).name; %pull out the 4 letter identifier

    [hdr, record] = edfread(files(i).name); %pull out data
    EEG1 = record(1,:);
    EEG2 = record(2,:);
    HPC = record(3,:)
    %EEG1 = filtfilt(b, a, record(1,:)); % Apply filter to raw EEG1
    %EEG2 = filtfilt(b, a, record(2,:)); % Apply filter to raw EEG1
    %HPC = filtfilt(b, a, record(3,:)); % Apply filter to raw EEG2

    clear record
    
    filepath='D:\JH\ephys\10.16_Combined_allCSTTrkBfiles_bydate\P32\'; % a subfolder inside the main fodler that holds the .edf (and produced) .mat files
    mousename=mouse_id{i};
    %save the contents to a .mat file
    save([filepath,files(i).name(1:12),'.mat'],'EEG1','EEG2', 'HPC', 'FS', 'mousename')
    clearvars -except files nfiles cutoffFreq filterOrder b a FS 
end


