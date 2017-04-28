#! /usr/bin/env python3


import argparse
from scipy.spatial.distance import pdist, squareform
from MyTools.points_io import point_reader
from matplotlib import pyplot as plt
import os
import numpy as np

__author__ = "Michał Kadlof <m.kadlof@cent.uw.edu.pl"


def draw(fname):
    points = point_reader(fname)
    n = len(points)
    base, _ = os.path.splitext(fname)
    ofname = base + '.pdf'
    k = 1000
    if len(points) > k:
        print('[WARN] lot of points - subsampling!')
        x = round(n/k)
        points = points[1::x]
    distances = squareform(pdist(points))
    plt.imshow(distances, cmap='viridis_r', origin='lower')
    if len(points) > k:
        plt.xticks(np.linspace(0, k, 6), np.linspace(0, n, 6).astype(np.int)+1)
        plt.yticks(np.linspace(0, k, 6), np.linspace(0, n, 6).astype(np.int)+1)
    plt.colorbar()
    plt.savefig(ofname)
    print('File {} saved.'.format(ofname))


def main():
    longer_help = """
    """

    parser = argparse.ArgumentParser(description="Draw distance map",
                                     epilog=longer_help, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("fname", help="PDB file name")
    args = parser.parse_args()
    draw(args.fname)


if __name__ == '__main__':
    main()