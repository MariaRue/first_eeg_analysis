
t=0.01:0.01:(60*1); % frame times
coh = randn(length(t),1); % vector of coherence frame by frame - replace with the real one

unf = coh; 
passbandfreq    = 10;
stopbandfreq    = 15;
passrip         = 1;
stopbandatten   = 60;
framerate       = 60; 
coherence_sd    = 0.05;
ft = designfilt('lowpassfir', 'PassbandFrequency', passbandfreq, 'StopbandFrequency', stopbandfreq,...
    'PassbandRipple', passrip, 'StopbandAttenuation', stopbandatten, 'SampleRate', framerate);
coh = filter(ft,coh); 



kernelDuration = 2000; % number of frames in sliding window

%kernel = ones(kernelDuration,1); %uniform kernel
kernel = fliplr(exp(-(1/2)*[1:100])); % exponential decay; halflife in frames is the number on the bottom of the fraction
%kernel = normpdf(0:kernelDuration, kernelDuration/2, 25); %Gaussian kernel

kernel = kernel./sum(kernel); %make area under curve be 1

coh_conv = conv(coh,kernel,'same');

figure; plot(t,coh,'k');
hold on; plot(t,coh_conv,'r');
hold off 

figure
[C,lag] = xcorr(coh_conv(1:100),coh_conv(1:100),'coeff'); 
plot(lag,C); 