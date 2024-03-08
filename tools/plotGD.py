import argparse
import gzip
from numpy import *
import matplotlib
import matplotlib.pyplot as plt
from scipy.special import logit
from scipy.special import expit
matplotlib.rcParams.update({'font.size': 10})

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot genotye Distribution data(s)")
    parser.add_argument("files", nargs='*', help="Data")
    parser.add_argument("--title", default="theta")
    parser.add_argument("--dMin", type=float, default=0.)
    parser.add_argument("--median",  action="store_true")
    parser.add_argument("--yMin",  type=float, default=1e-5)
    parser.add_argument("--yMax",  type=float, default=1e-3)
    parser.add_argument("--nPlots",  type=int, default=1)

    args = parser.parse_args()
    # get header and cols
    nFiles = len(args.files)

    ymi = [10., 10., 10., 10.]
    yma = [0., 0., 0., 0.]
    for i, file in enumerate(args.files):
        label = file.split(".txt.gz")[0].split("/")[-1]
        f         = gzip.open(file)
        idepths   = []
        imus      = []
        ithetas_g = []
        ithetas_r = []
        ihets     = []
        for j, key in enumerate(f.readline().split()):
            key = key.decode()
            if "depth" in key: idepths.append(j)
            if "mu" in key: imus.append(j)
            if "theta_g" in key: ithetas_g.append(j)
            if "theta_r" in key: ithetas_r.append(j)
            if key.endswith("het"): ihets.append(j)
        f.close()

        # get average and std
        data = genfromtxt(file, skip_header=1)
        if len(data.shape) == 1:
            data = array([data])

        depths = []
        for j in idepths:
            depths.append(data[:, j])

        mus = []
        for j in imus:
            mus.append(data[:, j])

        thetas_g = []
        for j in ithetas_g:
            thetas_g.append(data[:, j])

        thetas_r = []
        for j in ithetas_r:
            thetas_r.append(data[:, j])

        hets = []
        for j in ihets:
            hets.append(data[:, j])

        fmts= ["o-", "s-", "X-", "d-", "p-", "<-", "^-", ">-"]
        mks = [i for i in range(nFiles, 0, -1)]

        if args.median:
            print("using median")
            mdepths   = [nanmedian(d) for d in depths]
            mthetas_g = [nanmedian(t) for t in thetas_g]
            mthetas_r = [nanmedian(t) for t in thetas_r]
            mmus      = [nanmedian(t) for t in mus]
            mhets     = [nanmedian(t) for t in hets]
        else:
            print("using mean")
            mdepths   = [nanmean(d) for d in depths]
            mthetas_g = [nanmean(t) for t in thetas_g]
            mthetas_r = [nanmean(t) for t in thetas_r]
            mmus      = [nanmean(t) for t in mus]
            mhets     = [nanmean(t) for t in hets]

        sdepths   = [nanstd(d) for d in depths]
        sthetas_g = [nanstd(t) for t in thetas_g]
        sthetas_r = [nanstd(t) for t in thetas_r]
        smus      = [nanstd(t) for t in mus]
        shets     = [nanstd(t) for t in hets]

        print(label)
        print("depth:", mdepths)
        print("mu:", mmus)
        print("theta_r:", mthetas_r)
        print("theta_g:", mthetas_g)
        print("het:", mhets)



        ax1 = plt.subplot(411)
        plt.errorbar(mdepths, mthetas_g, xerr=sdepths, yerr=sthetas_g, fmt=fmts[i%len(fmts)], markersize=mks[i],linewidth=2, capsize=6, label=label)
        plt.xscale("log")
        plt.yscale("log")
        plt.ylabel(r"$\theta_g$")
        plt.ylim(args.yMin, args.yMax)
        ymi[0] = min(ymi[0], 10**(round(log10(min(mthetas_g))) - 1))
        yma[0] = max(yma[0], 10**(round(log10(max(mthetas_g))) + 1))
        plt.ylim(ymi[0], yma[0])
        plt.tick_params('x', labelbottom=False)

        plt.subplot(412, sharex=ax1)
        plt.errorbar(mdepths, mthetas_r, xerr=sdepths, yerr=sthetas_r, fmt=fmts[i%len(fmts)], markersize=mks[i],linewidth=2, capsize=6, label=label)
        plt.xscale("log")
        plt.yscale("log")
        plt.ylabel(r"$\theta_r$")
        ymi[1] = min(ymi[1], 10**(round(log10(min(mthetas_r))) - 1))
        yma[1] = max(yma[1], 10**(round(log10(max(mthetas_r))) + 1))
        plt.ylim(ymi[1], yma[1])
        plt.tick_params('x', labelbottom=False)

        plt.subplot(413, sharex=ax1)
        plt.errorbar(mdepths, mmus, xerr=sdepths, yerr=smus, fmt=fmts[i%len(fmts)], markersize=mks[i],linewidth=2, capsize=6, label=label)
        plt.xscale("log")
        plt.yscale("log")
        plt.ylabel(r"$\mu$")
        ymi[2] = min(ymi[2], 10**(round(log10(min(mmus))) - 1))
        yma[2] = max(yma[2], 10**(round(log10(max(mmus))) + 1))
        plt.ylim(ymi[2], yma[2])
        plt.tick_params('x', labelbottom=False)

        plt.subplot(414, sharex=ax1)
        plt.errorbar(mdepths, mhets, xerr=sdepths, yerr=shets, fmt=fmts[i%len(fmts)], markersize=mks[i],linewidth=2, capsize=6, label=label)
        plt.xscale("log")
        plt.yscale("log")
        plt.ylabel(r"heterozygosity")
        ymi[3] = min(ymi[3], 10**(round(log10(min(mhets))) - 1))
        yma[3] = max(yma[3], 10**(round(log10(max(mhets))) + 1))
        plt.ylim(ymi[3], yma[3])

        # All
        if i == nFiles - 1:
            plt.xlabel(r"Depth")
            plt.xlim(max(mdepths)*1.5, min(mdepths)/1.5)
            plt.xticks(mdepths, ["%2.2f"%(d) for d in mdepths])
            plt.legend()

    plt.show()
    #plt.savefig(args.files[0].split('_')[0] + ".png")
