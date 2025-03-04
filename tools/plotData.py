import argparse
import sys
import gzip
from numpy import *
import matplotlib
import matplotlib.pyplot as plt
from scipy.special import logit
from scipy.special import expit
matplotlib.rcParams.update({'font.size': 10})
tol_bright  = ['#4477AA', '#AA3377', '#228833', '#CCBB44', '#EE6677', '#66CCEE', '#BBBBBB', "#000000"]
tol_vibrant = ['#EE7733', '#0077BB', '#33BBEE', '#EE3377', '#CC3311', '#009988', '#BBBBBB']
IKRK        = ["#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]
basic       = ["#000000", "#FF0000", "#0000FF", "#FFA500", "#008000", "#808080", "#800080", "#008080"]

col  = tol_bright

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot columns of table")
    parser.add_argument("file", type=argparse.FileType('r'), nargs='?', default=sys.stdin, help="Transitions File")
    parser.add_argument("-x", type=int, default=1, help="x-column")
    parser.add_argument("-y", nargs='+', type=int, default=[2], help="y-column")
    parser.add_argument("--yMax", type=float, default=0., help="x-column")
    parser.add_argument("--logx", action="store_true", help="logx")
    parser.add_argument("--logy", action="store_true", help="logy")
    parser.add_argument("--vLine", "-v", type=int)
    parser.add_argument("--format", "-f", default="-")
    parser.add_argument("--grep", "-g", default="")

    args = parser.parse_args()
    header = args.file.readline().split()

    lines = filter(lambda line: args.grep in line, args.file)
    data = genfromtxt(lines)
    for y in args.y:
        plt.plot(data[:, args.x - 1], data[:, y - 1], args.format, label=header[y-1])

    if args.vLine is not None:
        plt.axvline(args.vLine, ls=":", color="k")

    if args.logx:
        plt.xscale("log")
    if args.logy:
        plt.yscale("log")


    if (args.yMax != 0.):
        plt.ylim(0, args.yMax)
    plt.legend()
    plt.xlabel(header[args.x-1])

    plt.show()
