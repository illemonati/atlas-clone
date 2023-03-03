import argparse
from numpy import *
import matplotlib
import matplotlib.pyplot as plt
from scipy.special import logit
from scipy.special import expit
matplotlib.rcParams.update({'font.size': 10})

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot Histogram")
    parser.add_argument("files", nargs='*', help="histogram-files")

    args = parser.parse_args()
    L = len(args.files)

    for fn in args.files:
        f = open(fn)
        header = f.readline().split() # skipfirst
        xname = header[1]
        yname = header[2]

        spl = f.readline().split()
        readgroup = spl[0]
        xs = [float(spl[1])]
        ys = [float(spl[2])]

        for l in f:
            spl = l.split()
            if spl[0] != readgroup:
                if L > 1:
                    plt.plot(xs, ys, label=fn)
                    break
                plt.plot(xs, ys, label=readgroup)
                readgroup = spl[0]
                xs = [float(spl[1])]
                ys = [float(spl[2])]
            xs.append(float(spl[1]))
            ys.append(float(spl[2]))

        if L == 1:
            plt.plot(xs, ys, label=readgroup)
            plt.xlabel(xname)
    plt.ylabel(yname)
    plt.legend()
    plt.show()
