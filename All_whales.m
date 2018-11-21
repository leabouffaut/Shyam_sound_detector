%% Shyam's Pitch tracking algorithm
clearvars    %MATLAB2016
close all
clc
addpath ../../../FONCTIONS
addpath ../../DATA_FORMAT
addpath matlab/detector
addpath matlab/util

% Extraction parameters
% ALL WHALES
annee = 2013;
station = 48;
jour = 148;
fact =  60/60;
bordure = 0 ; %(min)
duree =  fact*60.00 + bordure; %(min)
heure = 12.27;%12.30;%18.08;

% observation parameters
[x,fs,name] = cutfile_multiple_stations('../../',station, annee, jour, heure, duree);
x = x-mean(x); x = x/max(abs(x));
tx = (1:length(x))/fs;

% Spectrogram
taille_fft = 512;
overlap = 80; % \% overlapping
[stft,f,t,p] = spectrogram(x,hann(taille_fft),round((overlap/100)*taille_fft),taille_fft,fs);

% Removal of the frequencies below 10 Hz (non useful fere)
a =find(f>=12,1); 
b = find(f>=35,1);
p = p(a:b,:);
f = f(a:b);

fig = figure ; 
subplot 211
% p_affich = 10*log10(p);
% p_affich = p_affich+100; p_affich1 = (p_affich/max(max(p_affich)))*60;
% image(t,f,p_affich1);

imagesc(t,f,10*log10(p))
axis xy; axis tight; 
title('Original observation spectrogram')
% xlabel('Time (s)');
ylabel('Frequency (Hz)');
ylim([0 50])
grid on

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


hold on
for i = 1:length(my_tracks)
    plot(my_tracks{i}(1,:),my_tracks{i}(2,:),'k','Linewidth',2)
end

subplot 212
hold on
for i = 1:length(my_tracks)
    scatter(my_tracks{i}(1,:),my_tracks{i}(2,:),40,my_tracks{i}(3,:),'.')
end
grid on
ylabel('Pitch track (Hz)')
title('Pitch')
xlabel('Time (s)');
xlim([0 max(t)]); ylim([0 50])
