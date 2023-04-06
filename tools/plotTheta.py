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
    totalMax = 0   
    totalMin = 0
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

        real = means[0]

        probs = r_[probs]
        means = r_[means]
        stds  = r_[stds]

        print(file, probs, means, stds)

        fmts= ["o-", "s-", "x-", "d-", "<-"]
        mks = [10, 8, 6, 4, 2]
        plt.errorbar(probs, means, yerr=stds, fmt=fmts[i], markersize=8,linewidth=2, capsize=6, label=file)

        plt.xscale("log")
        plt.yscale("log")

        totalMax = max(totalMax, fileMax)
        totalMin = min(totalMin, fileMin)
        plt.hlines(real, 0, 10, "k", "dashed")

    plt.ylim(totalMin, totalMax+(totalMax*0.3))
    plt.legend()
    plt.xlabel(r"Downsampling probability")
    plt.ylabel(r"$\theta$")
    plt.title(args.title)
    plt.show()
