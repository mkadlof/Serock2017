#!/usr/bin/env python3

import re
import argparse
import numpy as np
import os
import sys
from matplotlib import pyplot as plt


def find_number_of_beads(fname):
	with open(fname) as f:
		for line in f:
			m = re.match('TER', line)
			if m:
				return int(line.split()[1])-1


def find_number_of_frames(fname):
	with open(fname, 'r') as f:
		for line in reversed(open(fname).readlines()):
			if "MODEL" in line:
				return int(line.split()[1])

 
def point_reader(fname):
    """Read the points from PDB file format

        Args:
            fname (string) filename of single chromatin model in pdb file format

        Returns:
            (list) List of three floats tuples representing points in euclidean R^3
    """
    atoms = [i.strip() for i in open(fname) if re.search('^(ATOM|HETATM)', i)]
    points = []
    for i in atoms:
        x = float(i[30:38])
        y = float(i[38:46])
        z = float(i[46:54])
        points.append((x, y, z))
    return np.array(points)


def read_trajectory(fname):	
	# find out how many beads is in the model
	base, _ = os.path.splitext(fname)
	npy_file = base + '.npy'
	try:
		trj = np.load(npy_file)
	except FileNotFoundError:
		no_beads = find_number_of_beads(fname)
		steps = find_number_of_frames(fname)
		trj = point_reader(fname).reshape((steps, no_beads, 3))
		np.save(npy_file, trj)
	return trj


def calculate_r_ends(trj):
	steps = trj.shape[0]
	n_beads = trj.shape[1]
	dystanse = []
	for i in range(steps):
		x1, y1, z1 = trj[i][0]
		x2, y2, z2 = trj[i][-1]
		dystanse.append(np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2))
	dystanse = np.array(dystanse)
	return dystanse


def calculate_rg(trj):
	steps = trj.shape[0]
	n_beads = trj.shape[1]
	rgs= []
	for i in range(steps):
		r0_x, r0_y, r0_z = np.mean(trj[i], axis=0)
		rg = 0
		for j in range(n_beads):
			x,y,z = trj[i][j]
			rg += (r0_x-x)**2+(r0_y-y)**2+(r0_z-z)**2
		rgs.append(rg)
	return np.array(rgs)
		
parser = argparse.ArgumentParser(description="opis")
parser.add_argument('infile', help="pdb trajectory")
args = parser.parse_args()
trj = read_trajectory(args.infile)
trj = trj[::30]		# subsumpling - weź co 10-tą komórkę trajektorii
dystanse = calculate_r_ends(trj)

rg = calculate_rg(trj)

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

l1 = ax1.plot(rg, lw=1, label="$R_g$")
l2 = ax2.plot(dystanse, lw=1, color='#fd8d3c', label="$R_k$")
fig.legend(l1+l2, ('$R_g$', '$R_k$') )
ax1.set_ylabel('promień bezwładności ($R_g$)')
ax2.set_ylabel('dystans między końcami ($R_k$)')
plt.grid()
plt.show()


