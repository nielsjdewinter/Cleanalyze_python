import numpy as np
import scipy.io as sio
from scipy.signal import find_peaks
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess
# from hampel import hampel

# Load the data
data = sio.loadmat('G730_sub.mat')
input_cleanalyze = data['G730_sub']
specimen = 'G730_sub'

# SN - change mm to um
input_cleanalyze[:, 0] = input_cleanalyze[:, 0] * 1000

# 1). REPLACE OUTLIERS WITH MEDIAN
temp, noutl = hampel(input_cleanalyze[:, 1], window_size=3)
out1 = np.column_stack((input_cleanalyze[:, 0], temp))
result1 = f'{np.sum(noutl)} outliers were identified and replaced with the local median'
print(result1)

# 2). AVERAGE DUPLICATE X-DATA
unique_x, indices, counts = np.unique(input_cleanalyze[:, 0], return_inverse=True, return_counts=True)
out2 = np.zeros((len(unique_x), 2))
out2[:, 0] = unique_x
out2[:, 1] = np.bincount(indices, weights=input_cleanalyze[:, 1]) / counts

# 3). MAKE DATA IN X EQUAL SPACED
Xinterpol = np.arange(0, np.floor(np.max(out2[:, 0])) + 0.5, 0.5)
interp_func = interp1d(out2[:, 0], out2[:, 1], kind='linear')
out3 = np.column_stack((Xinterpol, interp_func(Xinterpol)))

# 4). DETREND DATA -SN
nan_rows = np.isnan(out3).any(axis=1)
out3_clean = out3[~nan_rows]

t = out3_clean[:, 0]
p = np.polyfit(t, out3_clean[:, 1], 3)
f_y = np.polyval(p, t)
detr_data = out3_clean[:, 1] - f_y

# 5). GAUSSIAN SMOOTHING (GS) OF EQUAL SPACED DATA
out4GS = gaussian_filter1d(detr_data, sigma=20)

# 6). LOCATE LOCAL MAXIMA IN Y-DATA
out4GS_with_x = np.column_stack((out3_clean[:, 0], out4GS))
pksMaxima, _ = find_peaks(out4GS_with_x[:, 1], distance=1250, height=0.4)

# 7). LOESS filter -SN
x = out3[:, 0]
y = out3[:, 1]
y_smooth = lowess(y, x, frac=0.05, return_sorted=False)

# 8). Peak detection on top of LOESS - CODE BY NIELS
w = 1250
n = len(y_smooth)
window_size = 2 * w + 1
y_max = np.maximum.accumulate(y_smooth)
y_smooth_trimmed = y_smooth[w:n-w]
y_max_trimmed = y_max[w:n-w]
delta = y_max_trimmed - y_smooth_trimmed
i_max = np.where(delta <= 0)[0] + w

results = {
    'x': x[i_max],
    'i': i_max,
    'y_hat': y_smooth
}

# 9). PLOT THE RESULTS SN
titleFontSize = 20
subtitleFontSize = 20
axisLabelFontSize = 18
legendFontSize = 18
axisNumberFontSize = 18

plt.figure()
plt.subplot(2, 1, 1)
plt.suptitle(f'Specimen data: {specimen}', fontsize=subtitleFontSize)
plt.plot(out3_clean[:, 0], detr_data, '-', color='[0 0.4470 0.7410]', linewidth=1.5)
plt.plot(out3_clean[:, 0], out4GS, color='[1 0 0]', linewidth=2.5)
plt.plot(out3_clean[pksMaxima, 0], out4GS[pksMaxima], '.k', markersize=12)
plt.title('Detrended and Gaussian smoothed data with local maxima', fontsize=titleFontSize)
plt.legend(['Detrended data', 'Smoothed data using Gaussian', 'Local maxima'], loc='NE', fontsize=legendFontSize)
plt.xlabel('younger <- depth [µm] -> older', fontsize=axisLabelFontSize)
plt.ylabel('88Sr/43Ca [mmol/mol]', fontsize=axisLabelFontSize)
plt.ylim([-2, 8])
plt.xlim([0, 16000])
plt.grid(True)

plt.subplot(2, 1, 2)
plt.plot(out3[:, 0], out3[:, 1], '-', color='[0 0.4470 0.7410]', linewidth=1.5)
plt.plot(x, y_smooth, color='[1 0 0]', linewidth=2.5)
plt.plot(results['x'], y_smooth[results['i']], '.k', markersize=8, linewidth=12)
plt.title('Loess smoothed data with local maxima', fontsize=titleFontSize)
plt.legend(['Raw data', 'Smoothed data using loess', 'Local maxima'], loc='NE', fontsize=legendFontSize)
plt.xlabel('younger <- depth [µm] -> older', fontsize=axisLabelFontSize)
plt.ylabel('88Sr/43Ca [mmol/mol]', fontsize=axisLabelFontSize)
plt.ylim([0, 8])
plt.xlim([0, 16000])
plt.grid(True)

plt.show()
