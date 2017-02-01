%%JH 12/19/2016

%Like the line length script, data must be saved in the .mat format, with
%the data for each channel saved as a separate variable. 
load('8892MU.mat');

%% Z-score the data and make an index of when the data is above several SD

fs=2000;
dt=1/fs;
%maxtime=12600;
%time=dt:dt:maxtime;
%time=dt:dt:12600;
% HPC=HPC(1:25200000);

test_EEG= HPC;
 
%test_EEG = detrend_LFP(test_EEG');
test_EEG = test_EEG';

%test_normalized = bsxfun(@rdivide, bsxfun(@minus, test_EEG, baseline_power), baseline_std);
test_normalized = zscore(test_EEG);
test_normalizedabs = abs(test_normalized);

peak_idx = find(test_normalizedabs >= 7);

%% Finding the TOTAL number of events- WORKS 

x = diff(peak_idx)==1;%finds where the data is discontinuous, converts the data to binary
num_events = sum(x==0);%finds the places where the data is zero
events= cell(1,num_events);%creates a cell array to store all of the events


indends = find(x == 0); %lists the places where the data is a zero
ind = 1; 

for k = 1:num_events
   events{k} = peak_idx(ind:indends(k));
   ind = indends(k)+1;
end

events= events(cellfun('length', events) >=10);
cellsize=length(events);

for k=1:cellsize
    lengthevent = cellfun(@length, events);
end
 
%%
    average=mean(lengthevent);
     average=average*(1/fs);
     
 for k=1:length(events)
     first{k} = events{1,k}(1,1);
end


%%
epoched_data=cell(length(first));

for k=1:length(first)
    %length(first)
    epoched_times=((first(1,k)-10000):1:(first(1,k)+10000));
%     epoched_data{k}= test_EEG((epoched_times{1,k}));
end

%%
figure

for i=1:length(pot_seiz)
    subplot(5,4,i)
    plot(pot_seiz{1,i})
end
%% Length of each event- WORKS
 
 for k=1:length(peak_idx)
 EEGcolnum=peak_idx(k,1);
 seizetime(k)= time(EEGcolnum,1);
 end
 
