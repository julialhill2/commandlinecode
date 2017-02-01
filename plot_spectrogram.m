%Script to produce spectrogram figures for all files in file folder
%this script should be called after edftomat.m, which converts the raw .edf
%files to .mat files that are used by this script.

%First part uses a sample file to view the EEG and plot its spectrogram
%second part creates a big for loop to plot and save the spectrogram for
%all files in the PTZ folder

%Usage:
%To use this script, you have to to have a main directory (called mainpath
%below) that holds all of the .m files. Inside this directory, there should
%be a folder for the raw .edf files, which will be used to store the
%converted .mat files.

%Change lines 29-30 below to set the mainpath and the filepath for your PC.

%Madineh Sarvestani 2015 : continuation of discussion at CSHL Neural Data
%Science
%% Setup Matlab 
%clear the workspace, close all figures, clear the command prompt
clear all
clc
close all
warning off %turns off annoying warning messages


%% Set paths
mainpath='D:\JH\ephys\8.16_EEGEMG_NOR\Group2\OF\'; %update this with your path
%mainpath='C:\Users\Madineh\Dropbox\Julie\'; %update this with your path
filepath='D:\JH\ephys\8.16_EEGEMG_NOR\Group2\OF\'; % a subfolder inside the main fodler that holds the .edf (and produced) .mat files
%must have a backslash for the folder otherwise you will get error 'index
%exceeds matrix dimensions'

cd(mainpath);filepath
addpath(genpath(mainpath));
%% Preprocess: Load the Data
%we'll start with a sample file (file #1) and go through the steps. Then
%the last part of this script creates a big for loop to repeat these steps
%for all files

files=dir([filepath,'*.mat']); %find all the .edf files in folder
nfiles=length(files);
filenum=1;

load(files(filenum).name);
%% to look at a specific time for an event, seizure, whatever
start_time=1;
stop_time=2;

EEG1=EEG1((start_time*60*FS):(stop_time*60*FS));
EEG2=EEG2((start_time*60*FS):(stop_time*60*FS));
EMG=EMG((start_time*60*FS):(stop_time*60*FS));

%% Preprocess: Setup timevector for plotting
dt=1/FS;
maxtime=length(EEG1)/FS;
time=dt:dt:maxtime;

%can adjust by changing the length of the EEG file

%% View the data in short time windows, understand the range of the EEG values and exclude values
%that are outside of a reasonable range or the known extremems of your
%amplifier
%view the time-series
starttime=1; %in seconds
winlength=60; %size of full window
movelength=5; %how much to slide the window by on button press

mousename=files(filenum).name(end-7:end-4);
EEGslider(mousename,EEG,starttime,winlength,movelength,FS);
xlim([0 50])

%% View the data in spectrogram format: 
%we'll plot both the time-series and the spectrogram of all three channels

mousename=files(filenum).name(end-7:end-4);

%plot time-series
f1=figure('Position',[100 200 900 500]); %fix the figure position and give it a handle so we can save later
subplot(2,3,1); plot(time,EEG1);  ylim([-2000 2000]); 
title([num2str(mousename),',EEG1']);
xlabel('Time(s)'); axis tight;
subplot(2,3,2); plot(time,EEG2); ylim([-2000 2000]);
title([num2str(mousename),',EEG2']);
xlabel('Time (s)'); axis tight;
subplot(2,3,3); plot(time,EMG); ylim([-2000 2000]); 
title([num2str(mousename),',EMG']);
%xlabel('Time (s)'); axis tight;

%view the spectrogram of all signals in short (1 s) windows
%this is the trade-off between resolution in time and resolution in
%frequency: We want to have time resolution (running FFT on the whole
%signal turns the entire 30 minute signal into one frequency plot, so you
%can't tell which frequencies occured when). So we'll split the signal into
%short 1 s windows. This gives us better time resolution but worse
%frequency resolution (just because there are less data points to estimate
%the fft).

winlength=FS*1; 
%winlenth=FS*1=1 second window
noverlap=floor(0.5*winlength);  %overlap half the window 
%this means don't skip all the way to the next window, but slide a little
%bit. This is good for reducing edge artifacts that  MXCohen talked about.

%produce spectrogram using built-in matlab functions
%the outputs are: 
    %%spectrgrm is the forurier transform results (nxm matrix consisting of
    %%power in n frequency bins and m time bins)
    %%f is the frequencies
    %%t is the new time, calculated to include overlap amount
    %%p1 is the power in each of n frequency and m timebins. This is the
    %%output you want.

nfft=2^nextpow2(winlength);
[spectrgrm1,f,t,p1]=spectrogram(EEG1,winlength,noverlap,nfft,FS);
[spectrgrm2,f,t,p2]=spectrogram(EEG2,winlength,noverlap,nfft,FS);
[spectrgrm3,f,t,p3]=spectrogram(EMG,winlength,noverlap,nfft,FS);

%only look at 0-60 Hz, since we filtered everything else out
fmaxind=find(f>60,1,'first');
pplot1=p1(1:fmaxind,:);
pplot2=p2(1:fmaxind,:);
pplot3=p3(1:fmaxind,:);
fplot=f(1:fmaxind);

%plot the spectra
%EEG1
subplot(2,3,4); imagesc(10*log10(abs(pplot1))); 
axis xy; axis tight; colormap(jet); 
%setup time labels
numlabels=4; lengthtime=length(t); timebins=floor(lengthtime/numlabels);
set(gca,'xtick',1:timebins:length(t)); set(gca,'xticklabels',round(t(1:timebins:length(t))));
xlabel('Time (s)'); ylabel('Frequency (Hz)');
ylim([0 50]);
%now EEEG2
subplot(2,3,5); imagesc(10*log10(abs(pplot2))); 
axis xy; axis tight; colormap(jet); 
set(gca,'xtick',1:timebins:length(t)); set(gca,'xticklabels',round(t(1:timebins:length(t))));
xlabel('Time (s)'); ylabel('Frequency (Hz)'); ylim([0 50]);
%now HPC
subplot(2,3,6); imagesc(10*log10(abs(pplot3))); 
axis xy; axis tight; colormap(jet); 
set(gca,'xtick',1:timebins:length(t)); set(gca,'xticklabels',round(t(1:timebins:length(t))));
xlabel('Time (s)'); ylabel('Frequency (Hz)'); ylim([0 50]);

% now save this plot and name it with the mousename and group
% i always save matlab figures in two formats: fig, so I can open them in
% matlab and mess with labels, etc and jpg (or png, etc) so I can open them
% quickly in a folder

saveas(f1,['spectrogram ',mousename,'.fig'],'fig');
saveas(f1,['spectrogram ',mousename,'.jpg'],'jpg');

%% Now go through and produce the spectrogram for all 17 files using a big for loop
clear all 
close all
clc
warning off
addpath(genpath(cd));

files=dir('D:\JH\ephys\8.16_EEGEMG_NOR\Group2\OF\*.mat')
%files=dir('C:\Users\Madineh\Dropbox\Julie\PTZ\*.mat'); %find all the .edf files in folder
nfiles=length(files);

%setup parameters that won't change in the for loop
%setup time vector for plotting
FS=1000;
dt=1/FS;
%setup short window fourier transform parameters
winlength=FS*1; %1 second window
noverlap=floor(0.5*winlength);  %overlap half the window
nfft=2^nextpow2(winlength);

for i=2:nfiles
    display(['producing spectrogram plots for file ',num2str(i),' out of ', num2str(nfiles)]);

    load(files(i).name);
    
    maxtime=length(EEG1)/FS;
    time=dt:dt:maxtime;
    
    %run the spectrogram
    [~,f,t,p1]=spectrogram(EEG1,winlength,noverlap,nfft,FS);
    [~,f,t,p2]=spectrogram(EMG,winlength,noverlap,nfft,FS);
    [~,f,t,p3]=spectrogram(EEG2,winlength,noverlap,nfft,FS);
    
    %only look at 0-60 Hz, since we filtered everything else out
    fmaxind=find(f>60,1,'first');
    pplot1=p1(1:fmaxind,:);
    pplot2=p2(1:fmaxind,:);
    pplot3=p3(1:fmaxind,:);
    fplot=f(1:fmaxind);
    
    %plot the time domain
    f1=figure('Position',[100 200 900 500]); %fix the figure position and give it a handle so we can save later
    subplot(2,3,1); plot(time,EEG1);  ylim([-1000 1000]);
    title([num2str(mousename),',EEG1']);
    xlabel('Time(s)'); axis tight;
    subplot(2,3,2); plot(time,EMG); ylim([-1000 1000]);
    title([num2str(mousename),',EMG']);
    xlabel('Time (s)'); axis tight;
     subplot(2,3,3); plot(time,EEG2); ylim([-1000 1000]);
     title([num2str(mousename),',EEG2']);
     xlabel('Time (s)'); axis tight;
    
    %plot the spectra
    %EEG1
    subplot(2,3,4); imagesc(10*log10(abs(pplot1)));
    axis xy; axis tight; colormap(jet);
    %setup time labels
    numlabels=4; lengthtime=length(t); timebins=floor(lengthtime/numlabels);
    set(gca,'xtick',1:timebins:length(t)); set(gca,'xticklabel',round(t(1:timebins:length(t))));
    xlabel('Time (s)'); ylabel('Frequency (Hz)');
    %now EEG2
    subplot(2,3,5); imagesc(10*log10(abs(pplot2)));
    axis xy; axis tight; colormap(jet);
    set(gca,'xtick',1:timebins:length(t)); set(gca,'xticklabel',round(t(1:timebins:length(t))));
    xlabel('Time (s)'); ylabel('Frequency (Hz)');
     %now EMG
     subplot(2,3,6); imagesc(10*log10(abs(pplot3)));
     axis xy; axis tight; colormap(jet);
     set(gca,'xtick',1:timebins:length(t)); set(gca,'xticklabel',round(t(1:timebins:length(t))));
     xlabel('Time (s)'); ylabel('Frequency (Hz)');
     colorbar;
     
    %open a new results directory and plop figures in there
    %first make sure folder doesn't already exist
    if ~isempty([cd\'Results'])
    mkdir(cd,'Results');
    end
    % save
    savepath=[cd,'\Results\'];
    saveas(f1,[savepath,'spectrogram ',mousename,'.fig'],'fig');
    saveas(f1,[savepath,'spectrogram ',mousename,'.jpg'],'jpg');
    
    clearvars -except nfiles files i FS dt winlength noverlap nfft
    close all
    
end