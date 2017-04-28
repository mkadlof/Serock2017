# coding=utf-8

import csv
from matplotlib import pyplot as plt
import numpy as np

########################
#                      #
#   HELPER FUNCTIONS   #
#                      #
########################


def sizeof_fmt(num, suffix='B'):
    """Returns human readable size format (using binary prefixes)

        Args:
            num (float) : number to size
            suffix (str): a unit
    """
    for unit in ['', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi', 'Ei', 'Zi']:
        if abs(num) < 1024.0:
            return "{:3.1f} {}{}".format(num, unit, suffix)
        num /= 1024.0
    return "{:0.1f} {}{}".format(num, 'Yi', suffix)


def plot_data(state_fname, plot_fname):
    with open(state_fname) as csv_file:
        csv_reader = csv.DictReader(csv_file, fieldnames=['step', 'energy', 'temp'], delimiter=",")
        next(csv_reader, None)  # skip headers
        steps = []
        energies = []
        temp = []
        for row in csv_reader:
            steps.append(int(row['step']))
            energies.append(float(row['energy']))
            temp.append(float(row['temp']))
        steps = np.array(steps)
        energies = np.array(energies)
        temp = np.array(temp)
        N = 100
        temp_avg = np.convolve(temp, np.ones((N,))/N, mode='valid')
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()

        ax1.plot(steps, energies, lw=1, label="Potential energy")
        ax1.set_xlabel("Step")
        ax1.set_ylabel("Potential energy [kJ/mole]")

        ax2.plot(steps, temp, color='0.8', lw=0.5)
        ax2.plot(steps[:-N+1], temp_avg, color='red', lw=0.5, label='Temperature')
        ax2.set_ylabel("Temperature [K]")
        plt.savefig(plot_fname)
        print("Plot saved in  {} file".format(plot_fname))
