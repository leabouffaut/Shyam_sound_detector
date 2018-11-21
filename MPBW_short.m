%% Shyam's Pitch tracking algorithm
clearvars    %MATLAB2016
close all
clc
addpath ../../../FONCTIONS
addpath ../../DATA_FORMAT
addpath matlab/detector
addpath matlab/util

[x,fs,tx] = data_testbase_loading('MPBW');

% Spectrogram
taille_fft = 256;
overlap = 80; % \% overlapping
[stft,f,t,p] = spectrogram(x,hann(taille_fft),round((overlap/100)*taille_fft),taille_fft,fs);


%% Pitch detection

dT = (t(2)-t(1)); % Time delta for the spectrogram
dF = f(2)-f(1); %Frequency delta for the spectrogram
SNR_thresh = 4; % choose in the table in dt_SpectrogramRidgeDetector.m 
min_intensity = -Inf;
h_tracker = dt_SpectrogramRidgeTracker(dT, dF, f);
min_contour_len = 10; % mini nb of frames for a detection
max_contour_inactivity = 4; %compensate the Lloyds mirror effect

% % Easy way :
%  TestSpectrogramContourTracker(p(:,1:1000), dT, f, -Inf, 10, 4)

% Set a global variable called my_tracks where all the detection data are
% going to be put in the form of 1 listed item per track. Each track has 3
% rows: time, freq and power.
global my_tracks ;
my_tracks = []; 

% Set parameters for the class
% h_tracker.XXX = XXX function of the class described in
% dt_Spectrogram_Ridge_Detector

h_tracker.Set_threshold_value(SNR_thresh);
if ~isinf(min_intensity)
    h_tracker.Set_min_intensity(min_intensity);
end
h_tracker.SetMinContourLength(min_contour_len);
h_tracker.SetMaxContourInactivity(max_contour_inactivity);
h_tracker.SetTrackingStartTime(0);
h_tracker.SetCallback(@GatherTracksCB, 0);

% Analysis starts
h_tracker.ProcessFrames(p);
h_tracker.Flush();


%% PRINTING

fig = figure ; 
ax1 = subplot(2, 1, 1);
imagesc(t,f,10*log10(p)); % Plot original spectrogram
hold on
for i = 1:length(my_tracks) % Plot tracks
    plot(my_tracks{i}(1,:),my_tracks{i}(2,:),'k','Linewidth',2)
end
axis xy; axis tight; 
title('Original observation spectrogram')
ylabel('Frequency (Hz)');
ylim([0 50])
grid on


ax2 = subplot(2,1,2);
hold on
for i = 1:length(my_tracks)
    scatter(my_tracks{i}(1,:),my_tracks{i}(2,:),40,my_tracks{i}(3,:),'.')
end
grid on
ylabel('Pitch track (Hz)')
title('Pitch')
xlabel('Time (s)');
xlim([0 max(t)]); ylim([0 50])

o_c = get(ax1,'clim'); % color limit
set(ax2,'clim',o_c)
