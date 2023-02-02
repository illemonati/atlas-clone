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
    parser.add_argument("files", nargs='*', help="Data")
    parser.add_argument("--title", default="theta")

    args = parser.parse_args()
    real = -1
    sc   = 5

    # get header and cols
    for file in args.files:
        f     = gzip.open(file)
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
        data = genfromtxt(file, skip_header=1)
        means = []
        stds = []
        for i,c in enumerate(cols):
            means.append(mean(data[:, c]))
            stds.append(std(data[:, c]))

        if real < 0: real = means[0]

        probs = r_[probs]
        means = r_[means]
        stds  = r_[stds]

        plt.errorbar(probs, means, yerr=stds, fmt='o-', linewidth=2, capsize=6, label=file)

        plt.xscale("log")
        plt.yscale("log")
        sc = max(sc, max(means/real))
        sc = max(sc, max(real/means))

    sc = min(50, sc*2)
    plt.hlines(real, 0, 10, "k", "dashed")
    plt.ylim(real/sc, real*sc)
    plt.legend()
    plt.xlabel(r"Downsampling probability")
    plt.ylabel(r"$\theta$")
    plt.title(args.title)
    plt.show()
