#! /usr/bin/env python3

import sys
import argparse
import numpy as np
from numpy import pi, sin, cos, sqrt
from MyTools.points_io import save_points_as_pdb

__author__ = "Michał Kadlof <m.kadlof@cent.uw.edu.pl"


def line(n):
    points = []
    for i in range(n):
        points.append([i, 0, 0])
    return np.array(points)


def spiral_sphere(n, r, c):
    """based on: http://elib.mi.sanu.ac.rs/files/journals/vm/57/vmn57p2-10.pdf"""
    def calc_x(t): return sqrt(r**2-t**2)*cos(t/c)
    def calc_y(t): return sqrt(r**2-t**2)*sin(t/c)
    def calc_z(t): return t
    t_lin = np.linspace(-r, r, n)
    x = calc_x(t_lin).reshape((n, 1))
    y = calc_y(t_lin).reshape((n, 1))
    z = calc_z(t_lin).reshape((n, 1))
    points = np.concatenate((x, y, z), 1)
    return points


def baseball(n, r):
    """based on: http://paulbourke.net/geometry/baseball/"""
    b = pi/2
    a = 0.4     # shape parameter. Better not to change
    def calc_x(t): return r * sin(b - (b-a)*cos(t)) * cos(t/2. + a * sin(2*t))
    def calc_y(t): return r * sin(b - (b-a)*cos(t)) * sin(t/2. + a * sin(2*t))
    def calc_z(t): return r * cos(b - (b-a)*cos(t))
    t_lin = np.linspace(0, 4*pi-(4*pi)/n, n)
    x = calc_x(t_lin).reshape((n, 1))
    y = calc_y(t_lin).reshape((n, 1))
    z = calc_z(t_lin).reshape((n, 1))
    points = np.concatenate((x, y, z), 1)
    return points


def random_gas(n, d):
    x = np.random.uniform(-d, d, n).reshape((n, 1))
    y = np.random.uniform(-d, d, n).reshape((n, 1))
    z = np.random.uniform(-d, d, n).reshape((n, 1))
    points = np.concatenate((x, y, z), axis=1)
    return points


def random_walk(n):
    versors = np.random.uniform(-1, 1, 3*(n-1)).reshape(((n-1), 3))
    points = [[0, 0, 0]]
    for i in versors:
        points.append([points[-1][0]+i[0], points[-1][1]+i[1], points[-1][2] + i[2]])
    return np.array(points)


def build(args):
    if args.N > 9999:
        sys.exit("9999 beads is the limit! You requested for {}".format(args.N))
    error_msg = "Wrong number of parameters for {{}}. Please consult {0} -h".format(sys.argv[0])
    if args.type == "line":
        if len(args.params) != 0:
            sys.exit(error_msg.format(args.type))
        points = line(args.N)
        save_points_as_pdb(points, args.output)

    elif args.type == "spiral_sphere":
        if len(args.params) != 2:
            sys.exit(error_msg.format(args.type))
        r = args.params[0]
        c = args.params[1]
        points = spiral_sphere(args.N, r, c)
        save_points_as_pdb(points, args.output)

    elif args.type == "baseball":
        if len(args.params) != 1:
            sys.exit(error_msg.format(args.type))
        r = args.params[0]
        points = baseball(args.N, r)
        save_points_as_pdb(points, args.output)

    elif args.type == "random_gas":
        if len(args.params) != 1:
            sys.exit(error_msg.format(args.type))
        d = args.params[0]
        points = random_gas(args.N, d)
        save_points_as_pdb(points, args.output, render_connect=False)

    elif args.type == "random_walk":
        if len(args.params) != 0:
            sys.exit(error_msg.format(args.type))
        points = random_walk(args.N)
        save_points_as_pdb(points, args.output, render_connect=False)

    else:
        print("Wrong type {}. Please consult {} -h".format(args.type, sys.argv[0]))
        sys.exit(1)


def main():
    longer_help = """
    Additional params:
    line:
        no extra args
    spiral_sphere:
        r - radius of sphere (float)
        c - controls number of turns (float)
    baseball:
        r - radius of sphere (float)
    random_gas:
        d - dimension of box (-d, d)
    random_walk:
        no extra args
    """

    parser = argparse.ArgumentParser(description="Generate initial structure for polymer MD simulations",
                                     epilog=longer_help, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-c", '--center', type=bool, default=False, help="Center in PDB box? (5499, 5499, 4599)")
    parser.add_argument("-o", '--output', default="initial_structure.pdb", help="output PDB file name")
    parser.add_argument("type", help="line|spiral-sphere|baseball|random_gas")
    parser.add_argument("N", type=int, help="number of beads")
    parser.add_argument("params", metavar='params', type=float, nargs='*',
                        help="Additional params type dependent params")
    args = parser.parse_args()
    build(args)


if __name__ == '__main__':
    main()
