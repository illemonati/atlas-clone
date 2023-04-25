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
    parser.add_argument("--depth0", type=float, default=10.)

    args = parser.parse_args()
    # get header and cols
    for i, file in enumerate(args.files):
        f     = gzip.open(file)
        cols  = []
        heads = []
        probs =  []
        for j, key in enumerate(f.readline().split()):
            key = key.decode()
            if "theta_MLE" in key:
                cols.append(j)
                heads.append(key)
                if key[0] == "p":
                    probs.append(float(key[1:].split("_")[0]))
                else:
                    probs.append(1.)
        f.close()

        # get average and std
        data = genfromtxt(file, skip_header=1)
        if len(data.shape) == 1:
            data = array([data])
        means = []
        stds = []
        fileMax = 0
        fileMin = 0
        for c in cols:
            means.append(mean(data[:, c]))
            stds.append(std(data[:, c]))
            colMax = max(data[:, c])
            fileMax = max(fileMax, colMax)
            colMin = min(data[:, c])
            fileMin = min(fileMin, colMin)

        depths = r_[probs]*args.depth0
        means = r_[means]
        stds  = r_[stds]

        print(file, depths, means, stds)

        fmts= ["o-", "s-", "x-", "d-", "<-"]
        mks = [10, 8, 6, 4, 2]

        ax1 = plt.subplot(211)
        plt.errorbar(depths, means, yerr=stds, fmt=fmts[i], markersize=8,linewidth=2, capsize=6, label=file)
        plt.xscale("log")
        plt.yscale("log")
        plt.hlines(means[0], 0, 10, "k", "dashed")
        plt.subplot(212, sharex=ax1)
        plt.errorbar(depths, means/means[0], yerr=stds, fmt=fmts[i], markersize=8,linewidth=2, capsize=6, label=file)
        plt.xscale("log")
        plt.ylabel(r"$\theta/\theta_0$")
        plt.hlines(1, 0, 10, "k", "dashed")
        plt.xlabel(r"Depth")

    plt.subplot(211)
    plt.ylim(5e-4, 2e-3)
    plt.xlim(min(depths)/1.1, max(depths)*1.1)
    plt.legend()
    plt.title(args.title)
    plt.ylabel(r"$\theta$")
    plt.show()
    #plt.savefig(args.files[0].split('_')[0] + ".png")
