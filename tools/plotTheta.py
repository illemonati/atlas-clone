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
    parser.add_argument("--median",  action="store_true")
    parser.add_argument("--yMin",  type=float, default=1e-4)
    parser.add_argument("--yMax",  type=float, default=1e-2)
    parser.add_argument("--nPlots",  type=int, default=1)

    args = parser.parse_args()
    # get header and cols
    nFiles = len(args.files)
    ymi    = 10.
    yma    = 0.
    for i, file in enumerate(args.files):
        label = file.split("_theta.txt.gz")[0]
        f      = gzip.open(file)
        idepths = []
        ithetas = []
        il95s   = []
        iu95s   = []
        for j, key in enumerate(f.readline().split()):
            key = key.decode()
            if key.endswith("depth"): idepths.append(j)
            elif key.endswith("thetaMLE"): ithetas.append(j)
            elif key.endswith("thetaC95l"): il95s.append(j)
            elif key.endswith("thetaC95u"): iu95s.append(j)
        f.close()

        # get average and std
        data = genfromtxt(file, skip_header=1)
        if len(data.shape) == 1:
            data = array([data])

        depths = []
        for j in idepths:
            depths.append(data[:, j])

        thetas = []
        for j in ithetas:
            thetas.append(data[:, j])

        l95s = []
        for j in il95s:
            l95s.append(data[:, j])

        u95s = []
        for j in iu95s:
            u95s.append(data[:, j])


        fmts= ["o-", "s-", "X-", "d-", "p-", "<-", "^-", ">-"]
        mks = [10, 9, 8, 7, 6, 5, 4, 3, 2]

        for j in range(len(thetas)):
            iis = where(u95s[j] < 10.)
            thetas[j] = thetas[j][iis]
            print(j, len(l95s[j]), len(thetas[j]))

        if args.median:
            print("using median")
            mdepths = [nanmedian(d) for d in depths]
            mthetas = [nanmedian(t) for t in thetas]
        else:
            print("using mean")
            mdepths = [nanmean(d) for d in depths]
            mthetas = [nanmean(t) for t in thetas]

        sdepths = [nanstd(d) for d in depths]
        sthetas = [nanstd(t) for t in thetas]

        print(label)
        print("depth:", mdepths)
        print("theta:", mthetas)

        plt.errorbar(mdepths, mthetas, xerr=sdepths, yerr=sthetas, fmt=fmts[i], markersize=mks[i],linewidth=2, capsize=6, label=label)
        plt.hlines(mthetas[0], 0, 1.5*max(mdepths), plt.gca().lines[-1].get_color(), "dashed")
        print(plt.gca().lines[-1].get_color())
        plt.xscale("log")
        plt.yscale("log")
        plt.ylabel(r"$\theta$")
        plt.xlabel(r"Depth")
        plt.xlim(max(mdepths)*1.5, min(mdepths)/1.5)
        plt.xticks(mdepths, ["%2.2f"%(d) for d in mdepths])
        plt.legend()
        mthetas = r_[mthetas]
        mth_nz = mthetas[nonzero(mthetas)]
        ymi = min(ymi, 10**(round(log10(min(mth_nz))) - 1))
        yma = max(ymi, 10**(round(log10(max(mth_nz))) + 1))
        print(ymi, yma)
        plt.ylim(ymi, yma)

    plt.show()
    #plt.savefig(args.files[0].split('_')[0] + ".png")
