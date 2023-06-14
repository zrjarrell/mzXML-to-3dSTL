
from tkinter import filedialog
from pyteomics import mzxml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import surf2stl

class LCMSacquisition:
    def __init__(self, file_path):
        scanTimes = []
        spectrumFrames = []
        with mzxml.read(file_path) as reader:
            for spectrum in reader:
                scanTimes += [spectrum['retentionTime']]
                spectrumFrames += [pd.DataFrame(list(zip(spectrum['m/z array'], spectrum['intensity array'])), columns = ["mz", "intensity"])]
        scanTimes = np.array(scanTimes) * 60
        spectrumFrames = np.array(spectrumFrames, dtype=object)
        self.times = scanTimes
        self.spectra = spectrumFrames
    
    def roundTimes(self, sigFig=2):
        for time in range(0, len(self.times)):
            self.times[time] = round(10**sigFig * self.times[time]) / 10**sigFig

def getScan():
    file_path = filedialog.askopenfilename()
    newAcquisition = LCMSacquisition(file_path)
    return newAcquisition

def getSubset(analysis, mzWindow, timeWindow):
    subsetTimes = analysis.times[np.where((analysis.times > timeWindow[0]) & (analysis.times < timeWindow[1]))]
    subsetSpectra = analysis.spectra[np.where((analysis.times > timeWindow[0]) & (analysis.times < timeWindow[1]))]
    slicedSpectra = []
    for spectrum in subsetSpectra:
        slicedSpectra += [spectrum[(spectrum.loc[:,'mz'] > mzWindow[0]) & (spectrum.loc[:,'mz'] < mzWindow[1])].reset_index(drop=True)]
    return subsetTimes, slicedSpectra
    

###Inputs
mzWindow = [100, 101]
timeWindow = [50, 51]


qstd = getScan()
subsetTimes, subsetSpectra = getSubset(qstd,  mzWindow, timeWindow)



maxIntensity = 0
for spectrum in smallFrames:
    if max(spectrum.loc[:, 'intensity']) > maxIntensity:
        maxIntensity = max(spectrum.loc[:, 'intensity'])

printFrames = []
for spectrum in smallFrames:
    spectrum['relIntensity'] = spectrum['intensity'] / maxIntensity
    printFrames += [spectrum[spectrum['relIntensity'] > 0.0025]]

histoFrames = []
startBins = []
endBins = []
for i in range(0, 100):
    startBins += [146.5 + 0.01 * i]
    endBins += [146.5 + 0.01 * (i + 1) - 0.001]

for spectrum in printFrames:
    startIntensity = []
    endIntensity = []
    for i in range(0, 100):
        subset = spectrum[(spectrum['mz'] > startBins[i]) & (spectrum['mz'] < endBins[i])]
        if len(subset) > 0:
            startIntensity += [max(subset['relIntensity'])]
            endIntensity += [max(subset['relIntensity'])]
        else:
            startIntensity += [0.0025]
            endIntensity += [0.0025]
    startData = pd.DataFrame(list(zip(startBins, startIntensity)), columns = ["mz", "relIntensity"])
    endData = pd.DataFrame(list(zip(endBins, endIntensity)), columns = ["mz", "relIntensity"])
    histoFrames += [startData.append(endData).sort_values('mz', axis = 0).reset_index(drop = True)]


# testFrame = histoFrames[0]
# x = np.array(testFrame['mz'])
# y = np.linspace(0, 0.15, 2)
# X, Y = np.meshgrid(x, y)
# Z = np.array([testFrame['relIntensity'].to_list()] * 2)

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# ax.plot_surface(X, Y, Z, cmap=plt.cm.YlGnBu_r)

# plt.show()

# surf2stl.write('3d-1amu-2sec.stl', X, Y, Z)





# deltas = []
# for spectrum in printFrames:
#     masses = spectrum['mz'].to_list()
#     for mass in range(1, len(masses)):
#         deltas += [masses[mass]-masses[mass-1]]





x = np.array(histoFrames[0]['mz'])

Z = np.array([[0.0025]*200, [0.0025]*200, histoFrames[0]['relIntensity'].to_list(), histoFrames[0]['relIntensity'].to_list(), [0.0025]*200,
              [0.0025]*200, histoFrames[1]['relIntensity'].to_list(), histoFrames[1]['relIntensity'].to_list(), [0.0025]*200,
              [0.0025]*200, histoFrames[2]['relIntensity'].to_list(), histoFrames[2]['relIntensity'].to_list(), [0.0025]*200,
              [0.0025]*200, histoFrames[3]['relIntensity'].to_list(), histoFrames[3]['relIntensity'].to_list(), [0.0025]*200,
              [0.0025]*200, histoFrames[4]['relIntensity'].to_list(), histoFrames[4]['relIntensity'].to_list(), [0.0025]*200,
              [0.0025]*200, histoFrames[5]['relIntensity'].to_list(), histoFrames[5]['relIntensity'].to_list(), [0.0025]*200,
              [0.0025]*200, histoFrames[6]['relIntensity'].to_list(), histoFrames[6]['relIntensity'].to_list(), [0.0025]*200,
              [0.0025]*200, histoFrames[7]['relIntensity'].to_list(), histoFrames[7]['relIntensity'].to_list(), [0.0025]*200, [0.0025]*200])

y = np.array([round(subsetTimes[0]), subsetTimes[0] - 0.01, subsetTimes[0], subsetTimes[0] + 0.05, subsetTimes[0] + 0.06,
              subsetTimes[1] - 0.01, subsetTimes[1], subsetTimes[1] + 0.05, subsetTimes[1] + 0.06,
              subsetTimes[2] - 0.01, subsetTimes[2], subsetTimes[2] + 0.05, subsetTimes[2] + 0.06,
              subsetTimes[3] - 0.01, subsetTimes[3], subsetTimes[3] + 0.05, subsetTimes[3] + 0.06,
              subsetTimes[4] - 0.01, subsetTimes[4], subsetTimes[4] + 0.05, subsetTimes[4] + 0.06,
              subsetTimes[5] - 0.01, subsetTimes[5], subsetTimes[5] + 0.05, subsetTimes[5] + 0.06,
              subsetTimes[6] - 0.01, subsetTimes[6], subsetTimes[6] + 0.05, subsetTimes[6] + 0.06,
              subsetTimes[7] - 0.01, subsetTimes[7], subsetTimes[7] + 0.05, subsetTimes[7] + 0.06, round(subsetTimes[7])])
y = y / 2
#y = np.linspace(0, 1, len(Z))
X, Y = np.meshgrid(x, y)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(X, Y, Z, cmap=plt.cm.YlGnBu_r)

plt.show()





#whole map view
x = np.array(histoFrames[0]['mz'])

Z = np.array([[0.0025]*200, [0.0025]*200])
for i in histoFrames:
    Z = np.append(Z, [i['relIntensity']], axis=0)

y = np.linspace(0, 300, len(scanTimes)+2)
#y = np.linspace(0, 1, len(Z))
X, Y = np.meshgrid(x, y)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(X, Y, Z, cmap=plt.cm.YlGnBu_r)

plt.show()














subsetTimes = scanTimes[387:397]
subsetFrames = spectrumFrames[387:397]

for time in range(0, len(subsetTimes)):
    subsetTimes[time] = round(100 * subsetTimes[time]) / 100

smallFrames = []
for spectrum in subsetFrames:
    smallFrames += [spectrum[(spectrum.loc[:,'mz'] < 147.12) & (spectrum.loc[:, 'mz'] > 147.11)]]

maxIntensity = 0
for spectrum in smallFrames:
    if max(spectrum.loc[:, 'intensity']) > maxIntensity:
        maxIntensity = max(spectrum.loc[:, 'intensity'])

printFrames = []
for spectrum in smallFrames:
    spectrum['relIntensity'] = spectrum['intensity'] / maxIntensity
    printFrames += [spectrum[spectrum['relIntensity'] > 0.0025]]

histoFrames = []
startBins = []
endBins = []
for i in range(0, 100):
    startBins += [147.11 + 0.0001 * i]
    endBins += [147.11 + 0.0001 * (i + 1) - 0.00001]

for spectrum in printFrames:
    startIntensity = []
    endIntensity = []
    for i in range(0, 100):
        subset = spectrum[(spectrum['mz'] > startBins[i]) & (spectrum['mz'] < endBins[i])]
        if len(subset) > 0:
            startIntensity += [max(subset['relIntensity'])]
            endIntensity += [max(subset['relIntensity'])]
        else:
            startIntensity += [0.0025]
            endIntensity += [0.0025]
    startData = pd.DataFrame(list(zip(startBins, startIntensity)), columns = ["mz", "relIntensity"])
    endData = pd.DataFrame(list(zip(endBins, endIntensity)), columns = ["mz", "relIntensity"])
    histoFrames += [startData.append(endData).sort_values('mz', axis = 0).reset_index(drop = True)]



x = np.array(histoFrames[0]['mz'])

x = x * 100

Z = np.array([[0.0025]*200, [0.0025]*200, histoFrames[0]['relIntensity'].to_list(), histoFrames[0]['relIntensity'].to_list(), [0.0025]*200,
            [0.0025]*200])


y = np.array([round(subsetTimes[0]), subsetTimes[0] - 0.01, subsetTimes[0], subsetTimes[0] + 0.05, subsetTimes[0] + 0.06, 102.1])
#y = np.linspace(0, 1, len(Z))
X, Y = np.meshgrid(x, y)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(X, Y, Z, cmap=plt.cm.YlGnBu_r)

plt.show()


surf2stl.write('3d-1hundredth-amu-zoomed.stl', X, Y, Z)