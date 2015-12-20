# Read results into python of 2IS Kalman Filter/LDS and plot
import csv
import matplotlib.pyplot as plt

plt.figure()
with open("original-lisp-data.csv", 'r') as originalDataFile:
    with open("reconstructed-lisp-data.csv", 'r') as reconstructedDataFile:
        originalReader, reconstructedReader = csv.reader(originalDataFile), csv.reader(reconstructedDataFile)
        for i in range(len(originalReader)):
            plt.subplot(210 + i)
            nextOriginal
            plt.plot(range(len(originalReader[i])), originalReader[i], 'bo',
                     range(len(reconstructedReader[i])), reconstructedReader[i], 'r--')
    plt.show()