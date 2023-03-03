import argparse
from numpy import *
import matplotlib
import matplotlib.pyplot as plt
from scipy.special import logit
from scipy.special import expit
matplotlib.rcParams.update({'font.size': 10})

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot PMD Table")
    parser.add_argument("file", help="PMD-file")
    parser.add_argument("--cols", "-c", default="7 8 10 12", help="PMD-file")
    parser.add_argument("--pStart", "-p", default=1, type=int,help="PMD-file")

    args = parser.parse_args()
    f = open(args.file)
    lines = f.readlines()
    cols  = list(map(int, args.cols.split()))

    for c in cols:
        spl = lines[c].split()
        name = "_".join(spl[0:5])
        vs  = list(map(float, spl[5 + args.pStart-1:]))
        xs  = r_[0:len(vs)] + args.pStart - 1
        plt.plot(xs, vs, label=name)

    plt.xlabel(r"Position")
    plt.ylabel(r"Count")
    plt.legend()
    plt.show()
