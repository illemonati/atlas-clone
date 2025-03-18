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

    args = parser.parse_args()
    data = genfromtxt(args.file, skip_header=1)
    newData = []
    for l in data:
        lNew = l[1:]
        s = sum(lNew)
        if s > 0:
            lNew /= s
        newData.append(log(lNew))
    plt.imshow(newData)
    plt.colorbar()
    plt.show()
