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
    parser.add_argument("--dMin", type=float, default=0.)

    args = parser.parse_args()
    # get header and cols
    for i, file in enumerate(args.files):
        label = file.split("_theta.txt.gz")[0]
        f      = gzip.open(file)
        ithetas  = {} 
        idepths = {}
        for j, key in enumerate(f.readline().split()):
            key = key.decode()
            if "depth" in key:
                if key[0] == "p":
                    p = float(key[1:].split("_")[0])
                else:
                    p = 1.
                if not p in idepths:
                    idepths[p] = []
                idepths[p].append(j)
            if "theta_MLE" in key:
                if key[0] == "p":
                    p = float(key[1:].split("_")[0])
                else:
                    p = 1.
                if not p in ithetas:
                    ithetas[p] = []
                ithetas[p].append(j)
        f.close()

        # get average and std
        data = genfromtxt(file, skip_header=1)
        if len(data.shape) == 1:
            data = array([data])
        
        thetas  = []
        sthetas = []
        depths  = []
        sdepths = []
        for p in ithetas.keys():
            d = []
            for c in ithetas[p]:
                d = r_[d, data[:, c]]
            thetas.append(nanmean(d))
            sthetas.append(nanstd(d))
            d = []
            for c in idepths[p]:
                d = r_[d, data[:, c]]
            depths.append(nanmean(d))
            sdepths.append(nanstd(d))

        depths  = r_[depths]
        iis = where(depths[argsort(depths)[::-1]] > args.dMin)

        depths  = depths[iis]
        thetas  = r_[thetas][iis]
        sthetas = r_[sthetas][iis]
        sdepths = r_[sdepths][iis]



        print(file, ":");
        s = "depth = ["
        for j in range(len(depths)):
            s += "%4.4f +- %4.4f, "%(depths[j], sdepths[j])
        print(s[:-2] + "]")
        s = "theta = ["
        for j in range(len(thetas)):
            s += "%8.6f +- %8.6f, "%(thetas[j], sthetas[j])
        print(s[:-2] + "]")


        fmts= ["o-", "s-", "X-", "d-", "p-", "<-", "^-", ">-"]
        mks = [10, 9, 8, 7, 6, 5, 4, 3, 2]

        plt.errorbar(depths, thetas, xerr=sdepths, yerr=sthetas, fmt=fmts[i], markersize=8,linewidth=2, capsize=6, label=label)
        plt.xscale("log")
        plt.yscale("log")
        plt.hlines(thetas[0], 0, depths[0], "k", "dashed")
        plt.xlabel(r"Depth")
        plt.ylabel(r"$\theta/\theta_0$")

    #plt.ylim(5e-4, 2e-3)
    #plt.xlim(min(depths)/1.1, max(depths)*1.1)
    plt.legend()
    plt.title(args.title)
    plt.ylabel(r"$\theta$")
    plt.show()
    #plt.savefig(args.files[0].split('_')[0] + ".png")
