import argparse
import gzip
from numpy import *
import matplotlib
import matplotlib.pyplot as plt
from scipy.special import logit
from scipy.special import expit
matplotlib.rcParams.update({'font.size': 10})

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot theta data(s)")
    parser.add_argument("file", help="Data")

    args = parser.parse_args()

    # get header and cols
    f     = gzip.open(args.file)
    cols  = []
    heads = []
    probs =  []
    for i, key in enumerate(f.readline().split()):
        key = key.decode()
        if "theta_MLE" in key:
            cols.append(i)
            heads.append(key)
            probs.append(float(key[1:].split("_")[0]))
    f.close()

    # get average and std
    data = genfromtxt(args.file, skip_header=1)
    means = []
    stds = []
    for i,c in enumerate(cols):
        means.append(mean(data[:, c]))
        stds.append(std(data[:, c]))

    probs = r_[probs]
    means = r_[means]
    stds  = r_[stds]

    plt.errorbar(probs, means, yerr=stds, fmt='o', linewidth=2, capsize=6)

    real = means[0]

    ymin = min(means)
    ymax = max(means)
    dy   = max(stds)
    xmin = min(probs)
    xmax = max(probs)

    plt.hlines(real, 0, 10, "k", "dashed")

    plt.xscale("log")
    plt.xlim(xmin/2, xmax*2)

    plt.yscale("log")
    #plt.ylim(ymin - dy, ymax + dy)
    plt.ylim(ymin/2, ymax*2)

    plt.title(args.file)
    plt.xlabel(r"Downsampling probability")
    plt.ylabel(r"$\theta$")
    plt.show()
