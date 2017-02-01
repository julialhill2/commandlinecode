%% Line length script- written to analyze seizure data from CST/TrkB mutant mice

%Data must be in .mat format prior to use. File is loaded by double
%clicking on the file in the folder or by typing "load('xyz.mat')

%This file calculates line length based on Esteller 2001, subtracting
%successive data points, turning the result into absolute value, and then
%summing together all of the values. This is done within a 2 second window,
%and then the average line length for all of the bins for an animal is
%calculated as well as the St Dev. Then bins outside of this St Dev are
%identified. 

%1. Binning the specified data for analysis%%

fs=2000;
dt=1/fs;
total_t=length(EEG1);
timebin_length=fs*2;%(2 seconds);
num_bins=(total_t/timebin_length)-1; 
%number of bins to be created for the present data set-has to be 1 less
%than the total length 

bins=1:timebin_length:total_t;%creates a vector with the designated bins
bin_data=cell(1,num_bins);

test_EEG= EEG1;

%Loop that segments the data into bins based on the designated timebin
%length
for k=1:(length(bins)-1);
    a=bins(k);
    b=bins(k+1);
    datapoints=a:1:b;
    bin_data{k}=EEG1(datapoints);%this returns the data for each bin in a cell array
    
%     EEGstart=(EEG1(win));
%     EEGend=(EEG1(win+3999));
%     bin_data{k}=(EEGstart):(EEGend);
end

LLdist=length(datapoints);

%% Calculating line length- this part takes a few minutes

%allocating variables

summed_distances=cell(1,length(bin_data));
LL=zeros(1,length(bin_data));
final_sum=zeros(1,length(bin_data));

newdataval=zeros(1,4000);
firstval=zeros(1,4000);
nextval=zeros(1,4000);
distance_sum=zeros(1,4000);



%This loop does the line length measurement- finding the distance between
%each microvolt point in the binned data and containing those difference
%numbers in a cell array, where each cell contains the difference numbers
%for that bin of data. 

for k=1:length(bin_data);
    data=bin_data{k};
        
        for i=1:(length(data)-1)
            firstval(i)=data(1,i);
            nextval(i)=data(1,(i+1));
            diff=nextval(i)-firstval(i);
            diff=abs(diff);
            distance_sum(i)=diff;
        end
        
    summed_distances{k}=distance_sum;
        
        %summmed_distances{i}=summed_distances(i)
        %summed_distances{i}=(data(1,i)-(data(1,(i-1)));
    final_sum(k)=sum(summed_distances{k});
    LL(k)=(final_sum(k));

    %abs_bin_d{k}=abs(bin_data{k});%turns the EEG data into absolute values

    %sum_bin_d(k)=sum(abs_bin_d{k});%finds a sum of the absolute value for each timebin
    %LL(k)=(sum_bin_d(k))/(LLdist);%divides the summed absolute value by the #data points
    end

avg_LL=mean(LL);%finds the avg line length across all the bins in the data
std_LL=std(LL);%finds the standard deviation
twostd_LL=std_LL*2;%calculates a SD value for thresholding

thresh=avg_LL+twostd_LL;%creates a threshold 2 SD above the average line length
seiz=find(LL>thresh);%finds bins where the line length is 3 SD above thresh
noseiz=find(LL<thresh);%finds bins where the line length is below the thresh


%% Pulling out the data for potential seizures, bins where LL>thresh

EEGbins=bin_data;

for h=1:length(seiz); %this loop takes the bins with potential seizures and retrieves the data in cells
  seize_time=seiz(h);  
  pot_seiz{h}=EEGbins{seize_time}(1,:); %pot_seiz= cell array with the binned seizure data

end

for k=1:length(noseiz); %this loop takes the bins WITHOUT seizures and retrieves the data in cells (for power analysis)
    noseiz_time=noseiz(k);
    no_seiz{k}=EEGbins{noseiz_time}(1,:);
end

%% Plotting the seizures
%plots the first 20 instances above the threshold
tot=length(pot_seiz);

figure


for h=1:length(pot_seiz)
    subplot(5,4,h)
    plot(pot_seiz{1,h})
end


%% Spectra of potential seizures only
%Puts the data identified as above threshold into continuous data for power
%analysis 
totalseizdat=length(pot_seiz)*LLdist;
EEG_seizdat=zeros(1,totalseizdat);
EEG_seizdat=cell2mat(pot_seiz); %takes data from cell array and turns it into continous data for power spectra analysis 


window=FS*2;
noverlap=floor(0.5*window);
nfft=2^nextpow2(window);

[pxx,f] = pwelch(EEG_seizdat,window,noverlap,nfft,FS);
poten_seize_spectra = (10*log10(pxx));


%% Spectra of non-seizure data

total_no_seizdat=length(no_seiz)*LLdist;
EEG_noseiz_dat=zeros(1,total_no_seizdat);
EEG_noseiz_dat=cell2mat(no_seiz);

window=FS*2;
noverlap=floor(0.5*window);
nfft=2^nextpow2(window);

[pxx,f] = pwelch(EEG1,window,noverlap,nfft,FS);

normal_spectra = (10*log10(pxx));

%% Plot the normal data and the seizure data for the animal together
figure;
subplot(1,2,1)
plot(f,normal_spectra);
xlim([0 50]);
ylim([0 50]);
subplot(1,2,2)
plot(f,poten_seize_spectra);
xlim([0 50]);
ylim([0 50]);

