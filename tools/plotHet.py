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
    parser = argparse.ArgumentParser(description="Plot genotye Distribution data(s)")
    parser.add_argument("files", nargs='*', help="Data")
    parser.add_argument("--median", "-m", action="store_true")
    parser.add_argument("--relative",  "-r", action="store_true")
    parser.add_argument("--num",  "-n", type=int, default=1)
    parser.add_argument("--out",  "-o", default="")

    args     = parser.parse_args()
    nFiles   = len(args.files)
    nSamples = int(ceil(nFiles/args.num))

    if (nSamples) > len(col): sys.exit("Can only compare a maximum of %d samples!"%(len(col)))

    ymi = 10.
    yma = 0.
    xmi = 1000.
    xma = 0.

    for file in args.files:
        f = gzip.open(file)
        for key in f.readline().split():
            key = key.decode()
        f.close()

    for i, file in enumerate(args.files):
        label = file.split(".txt.gz")[0].split("/")[-1]
        f         = gzip.open(file)
        idepths   = []
        ihets     = []

        for j, key in enumerate(f.readline().split()):
            key = key.decode()
            if "depth" in key: idepths.append(j)
            if key.endswith("het"): ihets.append(j)
            if "expHetMLE" in key: ihets.append(j)
        f.close()

        # get average and std
        data = genfromtxt(file, skip_header=1)
        if len(data.shape) == 1:
            data = array([data])

        if args.relative:
            depth0 = data[0, idepths[0]]
            het0     = data[0, ihets[0]]
        else:
            depth0 = 1.
            het0     = 1.

        depths = []
        for j in idepths:
            depths.append(data[:, j]/depth0)

        hets = []
        for j in ihets:
            hets.append(data[:, j]/het0)

        if args.median:
            print("using median")
            mdepths   = r_[[nanmedian(d) for d in depths]]
            mhets     = r_[[nanmedian(h) for h in hets]]
        else:
            print("using mean")
            mdepths   = r_[[nanmean(d) for d in depths]]
            mhets     = r_[[nanmean(h) for h in hets]]

        sdepths   = r_[[nanstd(d) for d in depths]]
        shets     = r_[[nanstd(h) for h in hets]]
        
        print(label)
        print("depth:", mdepths)
        print("het:", mhets)

        xmi = min(min(mdepths), xmi)
        xma = max(max(mdepths), xma)

        fmts = ["o", "s", "X", "d", "p", "<", "^", ">"]
        lins = ["-", ":", "--"]
        mks  = [i for i in range(nFiles, 0, -1)]
        ax1  = plt.subplot(111)

        plt.errorbar(mdepths, mhets, color=col[i%nSamples], yerr=shets, fmt=fmts[i%nSamples] + lins[int(i/nSamples)], markersize=mks[i],linewidth=2, capsize=6, label=label)
        plt.hlines(mhets[0], 0, 1.5*max(mdepths), col[i%nSamples], "dashed")
        plt.xscale("log")
        plt.legend(ncols=2, borderaxespad=0.)
        if args.relative:
            plt.yscale("linear")
            yma = max(yma, max(mhets[nonzero(mhets)]))
            plt.ylim(0, 1.1*yma)
            plt.ylabel(r"het/het$_0$")

        else:
            plt.yscale("log")
            yma = max(yma, max(mhets[nonzero(mhets)]))
            ymi = min(ymi, min(mhets[nonzero(mhets)]))
            plt.ylim(min(ymi/2, yma/10), yma*2)
            plt.ylabel(r"het")

    # All
    plt.xlabel(r"Depth")
    xxs = r_[100,50,20,10,5,2,1,0.5,0.2,0.1,0.05,0.02,0.01,0.005,0.002,0.001,0.0005,0.0002,0.0001,0.00005,0.00002,0.00001]
    if xma > 1.1:
        plt.xticks(xxs, ["%2.2f"%(x) for x in xxs])
    else:
        plt.xticks(xxs, ["%.0e"%(x) for x in xxs])

    plt.xlim(xma*1.1, xmi/1.1)

    if args.out == "":
        plt.tight_layout()
        plt.show()
    else:
        fig = plt.gcf()
        fig.set_size_inches(9, 9)
        plt.tight_layout()
        plt.savefig(args.out, dpi=300)
