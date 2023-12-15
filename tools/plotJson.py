import argparse
import json
from numpy import *
import matplotlib
import matplotlib.pyplot as plt
from scipy.special import logit
from scipy.special import expit
matplotlib.rcParams.update({'font.size': 10})

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot theta data(s)")
    parser.add_argument("files", nargs='*', help="Data")
    parser.add_argument("--title", default="json")

    args = parser.parse_args()
    # get header and cols
    nFiles = len(args.files)
    iMax   = len(args.files) + 1
    for i, file in enumerate(args.files):
        print(file)
        fi = open(file)
        ji = json.load(fi)
        recal = ji[list(ji.keys())[0]]['recal']
        intercept = recal['intercept']
        print("intercept:", intercept)

        quality=recal['quality']
        quality = r_[quality[list(quality.keys())[0]]]
        quality[:,1] += intercept

        plt.subplot(311)
        plt.plot(quality[:,0], quality[:,1], linewidth=iMax/(i+1), label=file)
        plt.ylabel("logit")

        plt.subplot(312)
        plt.plot(quality[:,0], -10*log10(expit(quality[:,1])), linewidth=iMax/(i+1), label=file)
        plt.ylabel("quality")

        plt.subplot(313)
        plt.semilogy(quality[:,0], expit(quality[:,1]), linewidth=iMax/(i+1), label=file)
        plt.xlabel("Quality score")
        plt.ylabel("Probability")
    plt.subplot(311)
    plt.plot(quality[:,0], logit(10**(-quality[:,0]/10)), "k-", label="0")
    plt.subplot(312)
    plt.plot(quality[:,0], quality[:,0], "k-", label="0")
    plt.subplot(313)
    plt.plot(quality[:,0], 10**(-quality[:,0]/10), "k-", label="0")
    plt.legend()
    plt.show()
