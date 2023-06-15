
from pyteomics import mzxml
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D, axes3d
import surf2stl

class LCMSacquisition:
    def __init__(self, file_path):
        scanTimes = []
        spectrumFrames = []
        numSpectra = 0
        with mzxml.read(file_path) as reader:
            for spectrum in reader:
                if int(spectrum['id']) > numSpectra:
                    numSpectra = int(spectrum['id'])
        with mzxml.read(file_path) as reader:
            c = 1
            for spectrum in reader:
                print("Building spectra %d of %d." % (c, numSpectra), end="\r")
                scanTimes += [spectrum['retentionTime']]
                spectrumFrames += [pd.DataFrame(list(zip(spectrum['m/z array'], spectrum['intensity array'])), columns = ["mz", "intensity"])]
                c += 1
        scanTimes = np.array(scanTimes) * 60
        spectrumFrames = np.array(spectrumFrames, dtype=object)
        self.times = scanTimes
        self.spectra = spectrumFrames
    
    def roundTimes(self, sigFig=3):
        for time in range(0, len(self.times)):
            self.times[time] = properRound(self.times[time], sigFig)

def properRound(n, decimals=0):
    expoN = n * 10 ** decimals
    if abs(expoN) - abs(math.floor(expoN)) < 0.5:
        result =  math.floor(expoN) / 10 ** decimals
    else:
        result = math.ceil(expoN) / 10 ** decimals
    if decimals <= 0:
        return int(result)
    return result

def getSubset(analysis, mzWindow, timeWindow):
    subsetTimes = analysis.times[np.where((analysis.times > timeWindow[0]) & (analysis.times < timeWindow[1]))]
    subsetSpectra = analysis.spectra[np.where((analysis.times > timeWindow[0]) & (analysis.times < timeWindow[1]))]
    slicedSpectra = []
    for spectrum in subsetSpectra:
        slicedSpectra += [spectrum[(spectrum.loc[:,'mz'] > mzWindow[0]) & (spectrum.loc[:,'mz'] < mzWindow[1])].reset_index(drop=True)]
    return subsetTimes, slicedSpectra

def transformIntensity(subsetSpectra, method):
    if method == 'squareroot':
        for spectrum in subsetSpectra:
            transformed = np.sqrt(spectrum['intensity'])
            spectrum['intensity'] = transformed
    elif method == 'cuberoot':
        for spectrum in subsetSpectra:
            transformed = np.cbrt(spectrum['intensity'])
            spectrum['intensity'] = transformed
    elif method == 'log2':
        for spectrum in subsetSpectra:
            transformed = np.log2(spectrum['intensity'])
            spectrum['intensity'] = transformed
    elif method == 'log10':
        for spectrum in subsetSpectra:
            transformed = np.log10(spectrum['intensity'])
            spectrum['intensity'] = transformed
    else:
        print("Please provide a valid method: 'squareroot', 'cuberoot', 'log2', 'log10', or 'none'")

def getRelIntensity(subsetSpectra, relThreshold):
    maxIntensity = 0
    for spectrum in subsetSpectra:
        if max(spectrum.loc[:, 'intensity']) > maxIntensity:
            maxIntensity = max(spectrum.loc[:, 'intensity'])
    for spectrum in subsetSpectra:
        spectrum['relIntensity'] = spectrum['intensity'] / maxIntensity
        spectrum['relIntensity'][spectrum.loc[:,'relIntensity'] < relThreshold] = relThreshold
    return maxIntensity

def binMasses(subsetSpectra, numBins, mzWindow, relThreshold):
    binnedFrames = []
    bins = []
    stepSize = (mzWindow[1] - mzWindow[0]) / numBins
    for i in range(0, numBins):
        bins += [mzWindow[0] + stepSize * i, mzWindow[0] + stepSize * (i + 1) - stepSize / 1000]
    c = 1
    for spectrum in subsetSpectra:
        print("Binning spectrum %d of %d." % (c, len(subsetSpectra)), end = "\r")
        binnedIntensity = []
        for i in range(0, numBins):
            subset = spectrum[(spectrum['mz'] > bins[2*i]) & (spectrum['mz'] < bins[2*i+1])]
            if len(subset) > 0:
                binnedIntensity += [max(subset['relIntensity']), max(subset['relIntensity'])]
            else:
                binnedIntensity += [relThreshold, relThreshold]
        binnedData = pd.DataFrame(list(zip(bins, binnedIntensity)), columns = ["mz", "intensity"])
        binnedFrames += [binnedData]
        c += 1
    return binnedFrames

def plotSpectra(binnedFrames):
    x = np.array(binnedFrames[0]['mz'])
    y = np.array([subsetTimes[0]])
    Z = np.array([binnedFrames[0]['intensity'].to_list()])
    for i in range(1, len(binnedFrames)):
        Z = np.append(Z, [binnedFrames[i]['intensity'].to_list()], axis = 0)
        y = np.append(y, [subsetTimes[i]], axis = 0)
    X, Y = np.meshgrid(x, y)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z, cmap=plt.cm.YlGnBu_r)
    plt.show()

###Inputs
file_path = "/Users/zrj/Documents/research/scripts/qstd.mzXML"
mzWindow = [100, 200]
timeWindow = [50, 60]
relThreshold = 0.001
numBins = 2000
transformMethod = 'squareroot'

qstd = LCMSacquisition(file_path)

qstd.roundTimes()
subsetTimes, subsetSpectra = getSubset(qstd,  mzWindow, timeWindow)
if transformMethod != 'none':
    transformIntensity(subsetSpectra, transformMethod)

maxIntensity = getRelIntensity(subsetSpectra, relThreshold)
binnedFrames = binMasses(subsetSpectra, numBins, mzWindow, relThreshold)

plotSpectra(binnedFrames)

surf2stl.write('3d-10amu-10s.stl', X, Y, Z)






#old notes on generating one frame
""" x = np.array(histoFrames[0]['mz'])

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


surf2stl.write('3d-1hundredth-amu-zoomed.stl', X, Y, Z) """