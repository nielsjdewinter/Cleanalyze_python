import numpy as np
import scipy.io
import scipy.signal
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

def hampel_filter(input_series, window_size, n_sigmas=3):
    n = len(input_series)
    new_series = input_series.copy()
    k = 1.4826  # scale factor for Gaussian distribution
    indices = []

    for i in range((window_size), (n - window_size)):
        x0 = np.median(input_series[(i - window_size):(i + window_size)])
        S0 = k * np.median(np.abs(input_series[(i - window_size):(i + window_size)] - x0))
        if np.abs(input_series[i] - x0) > n_sigmas * S0:
            new_series[i] = x0
            indices.append(i)

    return new_series, indices

def gaussfilt_new(x, y, window_size):
    return gaussian_filter1d(y, window_size)

def DSQUARE(data, window_size):
    return np.gradient(np.gradient(data, window_size), window_size)

# Load data
mat = scipy.io.loadmat('G730_sub.mat')
input_cleanalyze = mat['G730_sub']
specimen = 'G730_sub'

# SN - change mm to um
input_cleanalyze[:, 0] = input_cleanalyze[:, 0] * 1000

# 1). REPLACE OUTLIERS WITH MEDIAN
temp, noutl = hampel_filter(input_cleanalyze[:, 1], 3)
out1 = np.column_stack((input_cleanalyze[:, 0], temp))
result1 = f'{len(noutl)} outliers were identified and replaced with the local median'
print(result1)

# 2). AVERAGE DUPLICATE X-DATA
unique_x, indices, counts = np.unique(input_cleanalyze[:, 0], return_inverse=True, return_counts=True)
out2 = np.column_stack((unique_x, np.bincount(indices, weights=input_cleanalyze[:, 1]) / counts))

# 3). MAKE DATA IN X EQUAL SPACED
Xinterpol = np.arange(0, np.floor(np.max(out2[:, 0])), 0.5)
temp = np.interp(Xinterpol, out2[:, 0], out2[:, 1])
out3 = np.column_stack((Xinterpol, temp))

# Replace NaN values with 0
out3[np.isnan(out3)] = 0

# 4). GAUSSIAN SMOOTHING (GS) OF EQUAL SPACED DATA
out4GS = np.column_stack((Xinterpol, gaussfilt_new(out3[:, 0], out3[:, 1], 20)))

# 5). LOCATE D-SQUARE MAXIMA
DSq = DSQUARE(out3[:, 1], 200)
pksDSq, _ = scipy.signal.find_peaks(DSq, distance=100, height=0.4)
locsDSq = pksDSq

# 6). CLUSTER REDUCTION D-SQUARE MAXIMA
clusterthreshold_microns = 300
clusterthreshold = clusterthreshold_microns / 0.5

numMaxima = len(locsDSq)
reducedMaxima = np.zeros(2 * numMaxima, dtype=int)
reducedMaximaIndex = 0

currentCluster = np.zeros(10, dtype=int)
currentClusterIndex = 0
currentCluster[currentClusterIndex] = locsDSq[0]

for i in range(1, len(locsDSq)):
    if locsDSq[i] - locsDSq[i - 1] < clusterthreshold:
        currentClusterIndex += 1
        currentCluster[currentClusterIndex] = locsDSq[i]
    else:
        reducedMaxima[reducedMaximaIndex:reducedMaximaIndex + 2] = [currentCluster[0], currentCluster[currentClusterIndex]]
        reducedMaximaIndex += 2
        currentClusterIndex = 0
        currentCluster[currentClusterIndex] = locsDSq[i]

if currentClusterIndex > 0:
    reducedMaxima[reducedMaximaIndex:reducedMaximaIndex + 2] = [currentCluster[0], currentCluster[currentClusterIndex]]

reducedMaxima = reducedMaxima[:reducedMaximaIndex]
reducedMaxima = np.unique(reducedMaxima)

# Define the font size
titleFontSize = 20
subtitleFontSize = 20
axisLabelFontSize = 18
legendFontSize = 18
axisNumberFontSize = 18
gMarkerSize = 6
gMarkerLineWidth = 2.5
roMarkerSize = 8
roMarkerLineWidth = 2

# 7). PLOT THE RESULTS
plt.figure()
plt.suptitle(f'Specimen data: {specimen}', fontsize=subtitleFontSize)

# Subplot 1: DSq maxima without data reduction
plt.subplot(2, 1, 1)
plt.plot(out3[:, 0], out3[:, 1], '-')
plt.plot(out4GS[:, 0], out4GS[:, 1], 'r-', linewidth=2)
plt.plot(out3[locsDSq, 0], out3[locsDSq, 1], 'g*', markersize=gMarkerSize, linewidth=gMarkerLineWidth)
plt.legend(['Original data', 'Smoothed data', 'D-Square maxima'], fontsize=legendFontSize)
plt.xlabel('younger <- distance [µm] -> older', fontsize=axisLabelFontSize)
plt.ylabel('88Sr/43Ca [mmol/mol]', fontsize=axisLabelFontSize)
plt.title('88Sr/43Ca data and D-Square maxima without reduction ws =, cutoff DSq =', fontsize=titleFontSize)
plt.ylim([0, 6])
plt.xlim([0, 16000])
plt.grid(True)

# Subplot 2: Reduced D-Square Maxima around Peaks
plt.subplot(2, 1, 2)
plt.plot(out3[:, 0], out3[:, 1], '-')
plt.plot(out4GS[:, 0], out4GS[:, 1], 'r-', linewidth=2)
plt.plot(out3[locsDSq, 0], out3[locsDSq, 1], 'g*', markersize=gMarkerSize, linewidth=gMarkerLineWidth)
plt.plot(out3[reducedMaxima, 0], out3[reducedMaxima, 1], 'ro', markersize=roMarkerSize, linewidth=roMarkerLineWidth)
plt.legend(['Original data', 'Smoothed data', 'All D-Square maxima', 'Reduced D-Square maxima'], fontsize=legendFontSize)
plt.xlabel('younger <- distance [µm] -> older', fontsize=axisLabelFontSize)
plt.ylabel('88Sr/43Ca ratio', fontsize=axisLabelFontSize)
plt.title('Reduced D-Square maxima', fontsize=titleFontSize)
plt.ylim([0, 8])
plt.xlim([0, 16000])
plt.grid(True)

plt.show()
