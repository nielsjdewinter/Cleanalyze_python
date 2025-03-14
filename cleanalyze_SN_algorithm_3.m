
% ==============================================================
% Dr. Frank Peeters
% Department of Earth Sciences
% Faculty of Science, Vrije Universiteit
% De Boelelaan 1085, 1081 HV Amsterdam
% The Netherlands. email: f.j.c.peeters@vu.nl
% =============================================================

% SN - Copy cleanalyze_20240220, changed parameters (see cleanalyze.m)
% Cleans and analyzes a two column array (x and y data) performing the
% following tasks:

% 1). Remove outliers using the hampel function
% 2). Average duplicate rows using the unique function
% 3). Convert x y data to equal spaced data in x
% 4). Gaussian smoothing
% 5). Identification of local maxima (above a given threshold)
% 6). Identification of D-Square maxima
% 7). Plot figure / results

clear, clc

wsOL = 3; % window size outlier detection (= windowsize hampel filter)

wsGS = 20; % window size for Gaussian smooting filter (1 sigma)

%wsLM = 2000; % window size for detecting local maxima
%peakthreshold = 2; % local maxima only identified above this value

wsDSq = 200; % window size for DSq analysis
cutoffDSq = 0.4; % cutoff level for D-Squeare filter

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

% Replace NaN values with 0
out3(isnan(out3)) = 0;
%% 4). GAUSSIAN SMOOTHING (GS) OF EQUAL SPACED DATA
out4GS = [Xinterpol,gaussfilt_new(out3(:,1),out3(:,2),wsGS)];

%% 5). LOCATE D-SQUARE MAXIMA
[DSq] = DSQUARE(out3(:,2),wsDSq);
[pksDSq,locsDSq] = findpeaks(DSq,'minpeakdistance',wsDSq/2,...
    'minpeakheight',cutoffDSq);

%% 6). CLUSTER REDUCTION D-SQUARE MAXIMA

% Define the clustering threshold in microns
clusterthreshold_microns = 300; % threshold converted to micrometers. 

% Convert the clustering threshold from microns to indices
clusterthreshold = clusterthreshold_microns / 0.5; 

% Initialize arrays to store the reduced D-Square maxima
numMaxima = length(locsDSq);
reducedMaxima = zeros(1, 2 * numMaxima);
reducedMaximaIndex = 1;

% Preallocate currentCluster with a reasonable maximum size
maxClusterSize = 10; % max size for storing DSq maxima indices in one single cluster before processing and resetting 
currentCluster = zeros(1, maxClusterSize);
currentClusterIndex = 1;
currentCluster(currentClusterIndex) = locsDSq(1);

for i = 2:length(locsDSq)
    if locsDSq(i) - locsDSq(i-1) < clusterthreshold
        currentClusterIndex = currentClusterIndex + 1;
        currentCluster(currentClusterIndex) = locsDSq(i);
    else
        % Process the current cluster
        reducedMaxima(reducedMaximaIndex:reducedMaximaIndex+1) = [currentCluster(1), currentCluster(currentClusterIndex)];
        reducedMaximaIndex = reducedMaximaIndex + 2;

        % Reset currentCluster for the next cluster
        currentClusterIndex = 1;
        currentCluster(currentClusterIndex) = locsDSq(i);
    end
end

% Add the last cluster
if currentClusterIndex > 0
    reducedMaxima(reducedMaximaIndex:reducedMaximaIndex+1) = [currentCluster(1), currentCluster(currentClusterIndex)];
end

% Remove unused preallocated space
reducedMaxima = reducedMaxima(1:reducedMaximaIndex);

% Remove duplicates and sort the reduced maxima
reducedMaxima = unique(reducedMaxima);

% Define the font size
titleFontSize = 20;
subtitleFontSize = 20;
axisLabelFontSize = 18;
legendFontSize = 18;
axisNumberFontSize = 18;
gMarkerSize = 6; % Define the marker size for g*
gMarkerLineWidth = 2.5; % Define the marker line width for g*
roMarkerSize = 8; % Define the marker size for ro
roMarkerLineWidth = 2; % Define the marker line width for ro

%% 7). PLOT THE RESULTS
figure;
sgtitle(['Specimen data:', specimen], 'Interpreter', 'none', 'FontSize', subtitleFontSize);

% Subplot 1: DSq maxima without data reduction
subplot(2,1,1);
plot(out3(:,1), out3(:,2), '-');
hold on;
plot(out4GS(:,1), out4GS(:,2), 'r-', 'LineWidth', 2);
plot(out3(locsDSq, 1), out3(locsDSq, 2), 'g*', 'MarkerSize', gMarkerSize, 'LineWidth', gMarkerLineWidth);
hold off;
legend('Original data', 'Smoothed data', 'D-Square maxima', 'FontSize', legendFontSize);
xlabel('younger <- distance [µm] -> older', 'FontSize', axisLabelFontSize);
ylabel('88Sr/43Ca [mmol/mol]', 'FontSize', axisLabelFontSize);
title('88Sr/43Ca data and D-Square maxima without reduction ws =, cutoff DSq =', 'FontSize', titleFontSize);
set(gca, 'YLim', [0 6], 'XLim', [0 16000], 'FontSize', axisNumberFontSize);
grid on;

% Subplot 2: Reduced D-Square Maxima around Peaks
subplot(2,1,2);
plot(out3(:,1), out3(:,2), '-'); % Plot original data
hold on;
plot(out4GS(:,1), out4GS(:,2), 'r-', 'LineWidth', 2); % Plot smoothed data
plot(out3(locsDSq, 1), out3(locsDSq, 2), 'g*', 'MarkerSize', gMarkerSize, 'LineWidth', gMarkerLineWidth); % Plot all D-Square maxima
plot(out3(reducedMaxima, 1), out3(reducedMaxima, 2), 'ro', 'MarkerEdgeColor', [1 0 0], 'MarkerSize', roMarkerSize, 'LineWidth', roMarkerLineWidth); % Plot reduced D-Square maxima
hold off;
legend('Original data', 'Smoothed data', 'All D-Square maxima', 'Reduced D-Square maxima', 'FontSize', legendFontSize);
xlabel('younger <- distance [µm] -> older', 'FontSize', axisLabelFontSize);
ylabel('88Sr/43Ca ratio', 'FontSize', axisLabelFontSize);
title('Reduced D-Square maxima', 'FontSize', titleFontSize);
set(gca, 'YLim', [0 8], 'XLim', [0 16000], 'FontSize', axisNumberFontSize);
grid on;
