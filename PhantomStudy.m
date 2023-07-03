% Generate Signals For Phantom Study

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  1: ON/OFF paradigm (Localiser)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate simple sinusoid at 20 Hz and alternate between 2s on 2s off for
% N trials. Beamform and use 4mm pseudo-T (maybe refine with 1mm??) to find
% position and orientation of phantom source dipole
close all
f_sample = 4000;
beta_amp = 1.2e-4; % amplitude in volts
sin_freq = 23;
on_dur = 2; %seconds
off_dur = 2;%on_dur;
trial_dur = on_dur + off_dur;
trial_time = linspace(0,trial_dur,f_sample*trial_dur);
N_trials = 100;

% win = hann(on_dur*f_sample)';
win = ones(1,on_dur*f_sample);

on_signal = beta_amp.*sin(2.*pi.*sin_freq.*trial_time(1:on_dur*f_sample)).*win;

off_signal = on_signal .*0;

trial_signal = [on_signal,off_signal];
trial_trig_signal = [ones(size(on_signal)),zeros(size(off_signal))];

%plot single trial signal, trigger and spectrum
figure
plot(trial_time,trial_signal)
hold on
plot(trial_time,trial_trig_signal)
figure
pwelch(trial_signal,[],[],[],f_sample)

outputSignal = repmat([trial_signal;trial_trig_signal],1,N_trials)';
all_time = linspace(0,length(outputSignal)./f_sample,length(outputSignal));
figure
plot(all_time,outputSignal)

% setup DAQ
dq = daq("ni");
dq.Rate = f_sample;
daq_name = "cDAQ2Mod1"; % depends on device name
addoutput(dq,daq_name , "ao0", "Voltage");
addoutput(dq, daq_name, "ao1", "Voltage");
write(dq, outputSignal)


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  2: faux "movie"/resting state paradigm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate signal with 1/f characteristic brain noise and pseudo
% oscillations in the alpha and beta range. Write this signal twice with a
% break for degaussing and recalibration etc.
close all

rng(42,'twister')

% A noise time-course is created by direct pole-placement. This creates a
% filter
% with an approximately 1/f frequency profile.
exp_length_seconds = 600;
f_sample = 4000;
time_vect = linspace(0,exp_length_seconds,exp_length_seconds*f_sample)';

% the value in poly is the filter root. it can vary between 0 < x < 1 where 0
% is white noise and 1 is an extremely sloped spectrum
% a is then the denominator polynomial of a digital filter which is used to
% create the signal

% relative amplitudes
alpha_amp = 155;
beta_amp = 80;
brain_noise_amp = 1;
norm_fac = sqrt(alpha_amp.^2 + beta_amp.^2 + brain_noise_amp.^2);

max_amp = 5e-4; %maximum amplitude in volts

%%%% brain noise
a = -poly(.99);
w_noise1 = randn(size(time_vect));
brain_noise = filter(1,a,w_noise1);
%%%% alpha signal
hp = 8;
lp = 12;
[b,a] = butter(3,2*[hp lp]/f_sample);
w_noise2 = randn(size(time_vect));
alpha_sig = [filtfilt(b,a,w_noise2)];


%%%% beta signal
hp = 13;
lp = 30;
[b,a] = butter(3,2*[hp lp]/f_sample);
w_noise3 = randn(size(time_vect));
beta_sig = [filtfilt(b,a,w_noise3)];


data = (brain_noise_amp/norm_fac).*brain_noise + ...
    (alpha_amp/norm_fac).*alpha_sig + ...
    (beta_amp/norm_fac).*beta_sig;

outputSignal = (max_amp.*data)./max(data);

daqreset
dq = daq("ni");
dq.Rate = f_sample;
addoutput(dq,daq_name, "ao0", "Voltage");
addoutput(dq, daq_name, "ao1", "Voltage");
% write(dq, [[0;outputSignal],[5;outputSignal*10000]])
% write(dq, [1e-3.*ones(320000,2)])

%
fftwin_s = 1;
[pxx,fxx] = pwelch(data,fftwin_s*f_sample,[],[],f_sample);

%
tlim = 5;

figure
tiledlayout('flow')
nexttile
plot(time_vect,brain_noise)
xlim([0,tlim])


nexttile
plot(time_vect,beta_sig)
xlim([0,tlim])

nexttile
plot(time_vect,outputSignal.*10000)
xlim([0,tlim])


nexttile
plot(fxx,pxx)
xlim([0,50])