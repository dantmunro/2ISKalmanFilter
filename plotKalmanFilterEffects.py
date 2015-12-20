# Read results into python of 2IS Kalman Filter/LDS and plot

# To obtain results, "kf" project must be first be run in Allegro CL
import csv
import matplotlib.pyplot as plt

NUM_PLOTS = 2

def isFloat(value):
    try:
        float(value)
        return True
    except:
        return False

plt.figure()
plt.suptitle("Comparisons of Noisy and Filtered Sine Waves")
with open("noisy-lisp-data.csv", 'r') as noisyDataFile:
    with open("filtered-lisp-data.csv", 'r') as filteredDataFile:
        noisyReader, filteredReader = csv.reader(noisyDataFile), csv.reader(filteredDataFile)
        for i in range(NUM_PLOTS):
            plt.subplot(211 + i)
            nextNoisy = [float(num) for num in noisyReader.next() if isFloat(num)]
            nextFiltered = [float(num) for num in filteredReader.next() if isFloat(num)]
            noisyPlot, = plt.plot(range(len(nextNoisy)), nextNoisy, 'bo', label="Noisy")
            filteredPlot, = plt.plot(range(len(nextFiltered)), nextFiltered, 'r--', label="Filtered")
            plt.legend([noisyPlot, filteredPlot], ["Noisy", "Filtered"])
    plt.show()
