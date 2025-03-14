
% ==============================================================
% Dr. Frank Peeters
% Department of Earth Sciences
% Faculty of Science, Vrije Universiteit
% De Boelelaan 1085, 1081 HV Amsterdam
% The Netherlands. email: f.j.c.peeters@vu.nl
% =============================================================

% SN - Copy cleanalyze_20240220, changed parameters (see cleanalyze.m)
% Original script developed by Dr. Frank Peeters
% Cleans and analyzes a two column array (x and y data) performing the
% following tasks:

% 1). Remove outliers using the hampel function
% 2). Average duplicate rows using the unique function
% 3). Convert x y data to equal spaced data in x
% 4). Detrend data - SN
% 5). Gaussian smoothing
% 6). Identification of local maxima (above a given threshold)
% 7). LOESS smoothing - SN
% 8). Peak detection on top of LOESS smoothing - N.J. de Winter
% 9). Plot figure / results - SN 

clear, clc

wsOL = 3; % window size outlier detection (= windowsize hampel filter)

wsGS = 20; % window size for Gaussian smooting filter (1 sigma)

wsLM = 2500; % window size for detecting local maxima
peakthreshold = 0.4; % local maxima only identified above this value

polynom = 3; % Order of polynomial for detrending the datasets

%wsDSq = 1000; % window size for DSq analysis
%cutoffDSq = 0.3; % cutoff level for D-Squeare filter

wsLO = 1250; % window size for peak detection in data smoothed with loess filter.
% half width window size. 

spLO = 0.05; % span loess filter, window size for loess filter, smoothing via local regression. 
% Span of 1% 


% load G732_sub.mat
% input_cleanalyze = G732_sub;
% specimen = 'G732_sub'; 

% load G737_sub.mat
% input_cleanalyze = G737_sub;
% specimen = 'G737_sub'; 

% load G717_sub.mat
% input_cleanalyze = G717_sub;
% specimen = 'G717_sub'; 

% load G781_sub.mat
% input_cleanalyze = G781_sub;
% specimen = 'G781_sub'; 

% load G753_sub.mat
% input_cleanalyze = G753_sub;
% specimen = 'G753_sub'; 

% load G767_sub.mat
% input_cleanalyze = G767_sub;
% specimen = 'G767_sub'; 

% load G664_sub.mat
% input_cleanalyze = G664_sub;
% specimen = 'G664_sub'; 

% load G770_sub.mat
% input_cleanalyze = G770_sub;
% specimen = 'G770_sub'; 

load G730_sub.mat
input_cleanalyze = G730_sub;
specimen = 'G730_sub'; 

% load G723_sub.mat
% input_cleanalyze = G723_sub;
% specimen = 'G723_sub'; 

% load G668_sub.mat
% input_cleanalyze = G668_sub;
% specimen = 'G668_sub'; 

% load G652_sub.mat
% input_cleanalyze = G652_sub;
% specimen = 'G652_sub';

% load G691_int.mat
% input_cleanalyze = G691_int;
% specimen = 'G691_int';

% load G728_int.mat
% input_cleanalyze = G728_int;
% specimen = 'G728_int';

% load G719_int.mat
% input_cleanalyze = G719_int;
% specimen = 'G719_int';

% load G779_int.mat
% input_cleanalyze = G779_int;
% specimen = 'G779_int';

% load B235_int.mat
% input_cleanalyze = B235_int;
% specimen = 'B235_int';

% load B232_int.mat
% input_cleanalyze = B232_int;
% specimen = 'B232_int';

% load input_cleanalyze_kokkel.mat
% input_cleanalyze = input_cleanalyze_kokkel;



% SN - change mm to um
%% 0). 
input_cleanalyze(:,1) = input_cleanalyze(:,1) * 1000;


%% 1). REPLACE OUTLIERS WITH MEDIAN
[temp,noutl] = hampel(input_cleanalyze(:,2),wsOL);
out1 = [input_cleanalyze(:,1), temp];
clear temp
result1 = [num2str(sum(noutl)),...
    ' outliers were identified and replaced with the local median'];
disp(result1)
% difference = input_cleanalyze(:,2)-out1(:,2);
% plot(difference,'+')

%% 2). AVERAGE DUPLICATE X-DATA
[out2,~,z] = unique(input_cleanalyze(:,1),'stable');
ncol = size(input_cleanalyze,2);
out2(:,2) = accumarray(z,input_cleanalyze(:,2),[],@mean);

%% 3). MAKE DATA IN X EQUAL SPACED
Xinterpol = 0:0.5:floor(max(out2(:,1)));
Xinterpol = Xinterpol';
temp = interp1(out2(:,1),out2(:,2),Xinterpol,'linear');
out3 = [Xinterpol, temp];
clear temp

%% 4). DETREND DATA -SN
% detrending of data with non-linear trend

% Find rows with NaN values
nan_rows = any(isnan(out3), 2);

% Remove rows with NaN values
out3_clean = out3(~nan_rows, :);
 
% Detrend using polynomial fitting
polynom; % Order of polynomial
t = out3_clean(:,1); % Use x-values of out3_clean as time vector
[p, ~, mu] = polyfit(t, out3_clean(:,2), polynom); % Fit an x degree polynomial
f_y = polyval(p, t, [], mu); % Evaluate the polynomial

% Detrended data
detr_data = out3_clean(:,2) - f_y;

%% 5). GAUSSIAN SMOOTHING (GS) OF EQUAL SPACED DATA
% Gaussian smoothing of data without detrending - original
% out4GS = [Xinterpol,gaussfilt_new(out3_clean(:,1),detr_data(:,1),wsGS)];

% Gaussian Smoothing of Detrended Data (SN)
out4GS = gaussfilt_new(out3_clean(:,1), detr_data, wsGS);

%% 6). LOCATE LOCAL MAXIMA IN Y-DATA 
% Adjusted SN
out4GS_with_x = [out3_clean(:, 1), out4GS]; % make out4GS 29851 x 2 instead of x 1. Adding the x-values from out3_clean.

% Find local maxima in detrended and smoothed data
[pksMaxima,locsMaxima] = findpeaks(out4GS_with_x(:,2),'minpeakdistance', wsLM/2,...
    'minpeakheight',peakthreshold); 

%% 7). LOESS filter -SN

% identifying x, y, w for peak detection, and span for loess smoothing
% (loess = smoothing via local regression)

x = out3(:,1);
y = out3(:,2);
wsLO;
span = spLO;

% Smoothing data using loess smoothing with a variable span
y_smooth = smooth(x, y, spLO, 'loess');

%% 8). Peak detection on top of LOESS - CODE BY NIELS
% Line-by-line translation of the R function made based on the blogpost
% Give the half-width of the window (w)
% x = depth; y = Sr/Ca;
% depth = max(out1);
% y = max(out2);

x = out3(:,1); % depth
y = out3(:,2); % Sr/Ca
 
w = wsLO; % Adjust to the right half-width.

% % Find the total length of the dataset
n = length(y_smooth);

% This calculates the full window size
window_size = 2 * w + 1;

% Apply the moving maximum function with a centered window
y_max = movmax(y_smooth, [w w]);

% Remove the first w and last w elements from y_smooth and y_max
y_smooth_trimmed = y_smooth((w+1):(n-w));
y_max_trimmed = y_max((w+1):(n-w));
 
% Calculate delta (the difference between y_smooth and y_max)
delta = y_max_trimmed - y_smooth_trimmed;

% Find indices where delta is less than or equal to zero
i_max = find(delta <= 0) + w;

% Create a structure to store the results
results.x = x(i_max);
results.i = i_max;
results.y_hat = y_smooth;
 
% The "results" array contains the location of the peaks (x), the index of the peaks (i) and the smoothed y value belonging to the peaks

%% 9). PLOT THE RESULTS SN

% Plot detrended and Gaussian smoothed data with local maxima
figure;
subplot(2,1,1);
sgtitle(['Specimen data:', specimen], 'Interpreter', 'none', 'FontSize', subtitleFontSize);
plot(out3_clean(:,1), detr_data, '-', 'Color', '[0 0.4470 0.7410]', 'LineWidth', 1.5); % Detrended data
hold on;
plot(out3_clean(:,1), out4GS, 'Color', '[1 0 0]', 'LineWidth', 2.5); % Smoothed data using Gaussian filter
plot(out3_clean(locsMaxima, 1), out4GS(locsMaxima), '.k', 'MarkerSize', 12); % Local maxima
hold off;
title('Detrended and Gaussian smoothed data with local maxima wsGS =, wsLM =, peak cutoff =', 'FontSize', titleFontSize);
legend('Detrended data', 'Smoothed data using ''Gaussian''', 'Local maxima', 'Location', 'NE', 'FontSize', legendFontSize);
xlabel('younger <- depth [µm] -> older', 'FontSize', axisLabelFontSize);
ylabel('88Sr/43Ca [mmol/mol]', 'FontSize', axisLabelFontSize);
set(gca, 'YLim', [-2 8], 'XLim', [0 16000], 'FontSize', axisNumberFontSize);
grid on;

% Plotting the raw data, smoothed LOESS data, and the detected peaks
figure;
subplot(2,1,2);
sgtitle(['Specimen data:', specimen], 'Interpreter', 'none', 'FontSize', subtitleFontSize);
plot(out3(:,1), out3(:,2), '-', 'Color','[0 0.4470 0.7410]', 'LineWidth', 1.5);
hold on;
plot(x, y_smooth, 'Color', '[1 0 0]', 'LineWidth', 2.5); % plotting the smoothed data
plot(results.x, y_smooth(results.i), '.k', 'MarkerSize', 8, 'LineWidth', 12); % Plotting the peaks
hold off;
title('Loess smoothed data with local maxima wsLS =, wsLM = ', 'FontSize', titleFontSize);
legend('Raw data', 'Smoothed data using ''loess''', 'Local maxima', 'Location', 'NE', 'FontSize', legendFontSize);
xlabel('younger <- depth [µm] -> older', 'FontSize', axisLabelFontSize);
ylabel('88Sr/43Ca [mmol/mol]', 'FontSize', axisLabelFontSize);
set(gca, 'YLim', [0 8], 'XLim', [0 16000], 'FontSize', axisNumberFontSize);
grid on;

% Define the font size
titleFontSize = 20;
subtitleFontSize = 20;
axisLabelFontSize = 18;
legendFontSize = 18;
axisNumberFontSize = 18;
