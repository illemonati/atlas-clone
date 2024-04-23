import argparse
import gzip
from numpy import *
import matplotlib
import matplotlib.pyplot as plt
from scipy.special import logit
from scipy.special import expit
matplotlib.rcParams.update({'font.size': 10})
tol_bright  = ['#4477AA', '#AA3377', '#228833', '#CCBB44', '#EE6677', '#66CCEE', '#BBBBBB']
tol_vibrant = ['#EE7733', '#0077BB', '#33BBEE', '#EE3377', '#CC3311', '#009988', '#BBBBBB']
IKRK        = ["#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]
basic       = ["#000000", "#FF0000", "#0000FF", "#FFA500", "#008000", "#808080", "#800080", "#008080"]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot genotye Distribution data(s)")
    parser.add_argument("files", nargs='*', help="Data")
    parser.add_argument("--median", "-m", action="store_true")
    parser.add_argument("--relative",  "-r", action="store_true")
    parser.add_argument("--one",  "-o", action="store_true")

    args = parser.parse_args()
    # get header and cols
    nFiles = len(args.files)

    ymi = [10., 10., 10., 10.]
    yma = [0., 0., 0., 0.]
    xmi = 1000.
    xma = 0.

    hky85 = False #at least one file is hky85
    if not args.one:
        for file in args.files:
            f = gzip.open(file)
            for key in f.readline().split():
                key = key.decode()
                if "mu" in key: hky85 = True
            f.close()

    for i, file in enumerate(args.files):
        label = file.split(".txt.gz")[0].split("/")[-1]
        f         = gzip.open(file)
        idepths   = []
        imus      = []
        ithetas_g = []
        ithetas_r = []
        ithetas_f = []

        for j, key in enumerate(f.readline().split()):
            key = key.decode()
            if "depth" in key: idepths.append(j)
            if "mu" in key: imus.append(j)
            if "theta_g" in key: ithetas_g.append(j)
            if "theta_r" in key: ithetas_r.append(j)
            if "thetaMLE" in key: ithetas_f.append(j)
        f.close()

        # get average and std
        data = genfromtxt(file, skip_header=1)
        if len(data.shape) == 1:
            data = array([data])

        depths = []
        if args.relative:
            depth0 = data[0, idepths[0]]
        else:
            depth0 = 1.
        for j in idepths:
            depths.append(data[:, j]/depth0)

        hky85_i = len(imus) > 0

        if hky85_i:
            label += "_HKY85"
            if args.relative:
                mu0      = data[0,imus[0]]
                theta0_g = data[0,ithetas_g[0]]
                theta0_r = data[0,ithetas_r[0]]
            else:
                mu0      = 1.
                theta0_g = 1.
                theta0_r = 1.

            mus = []
            for j in imus:
                mus.append(data[:, j]/mu0)

            thetas_g = []
            for j in ithetas_g:
                thetas_g.append(data[:, j]/theta0_g)

            thetas_r = []
            for j in ithetas_r:
                thetas_r.append(data[:, j]/theta0_r)
        else:
            label += "_Fels"
            if args.relative:
                theta0 = data[0,ithetas_f[0]]
            else:
                theta0 = 1.

            thetas_g = []
            for j in ithetas_f:
                thetas_g.append(data[:, j]/theta0)

        fmts = ["o-", "s-", "X-", "d-", "p-", "<-", "^-", ">-"]
        mks  = [i for i in range(nFiles, 0, -1)]
        col  = tol_bright

        if args.median:
            print("using median")
            mdepths   = r_[[nanmedian(d) for d in depths]]
            mthetas_g = r_[[nanmedian(t) for t in thetas_g]]
            if hky85_i:
                mthetas_r = r_[[nanmedian(t) for t in thetas_r]]
                mmus      = r_[[nanmedian(t) for t in mus]]
        else:
            print("using mean")
            mdepths   = r_[[nanmean(d) for d in depths]]
            mthetas_g = r_[[nanmean(t) for t in thetas_g]]
            if hky85_i:
                mthetas_r = r_[[nanmean(t) for t in thetas_r]]
                mmus      = r_[[nanmean(t) for t in mus]]

        sdepths   = r_[[nanstd(d) for d in depths]]
        sthetas_g = r_[[nanstd(t) for t in thetas_g]]
        if hky85_i:
            sthetas_r = r_[[nanstd(t) for t in thetas_r]]
            smus      = r_[[nanstd(t) for t in mus]]
        
        if hky85_i:
            print(label)
            print("depth:", mdepths)
            print("theta_g:", mthetas_g)
            print("theta_r:", mthetas_r)
            print("mu:", mmus)
        else:
            print(label)
            print("depth:", mdepths)
            print("theta_f:", mthetas_g)

        xmi = min(min(mdepths), xmi)
        xma = max(max(mdepths), xma)

        if hky85: ax1 = plt.subplot(311)
        else:     ax1 = plt.subplot(111)
        l = plt.errorbar(mdepths, mthetas_g, color=col[i], xerr=sdepths, yerr=sthetas_g, fmt=fmts[i%len(fmts)], markersize=mks[i],linewidth=2, capsize=6, label=label)
        plt.hlines(mthetas_g[0], 0, 1.5*max(mdepths), col[i], "dashed")
        plt.xscale("log")
        plt.legend(ncols=2, borderaxespad=0.)
        if args.relative:
            plt.yscale("linear")
            plt.ylim(0, 1.5)
            if hky85: plt.ylabel(r"$\theta_{f/g}/\theta_0$")
            else: plt.ylabel(r"$\theta_f/\theta_0$")
        else:
            plt.yscale("log")
            mas    = mthetas_g + sthetas_g
            mis    = mthetas_g - sthetas_g
            yma[0] = max(yma[0], max(mas[nonzero(mas)])*1.1)
            ymi[0] = max(yma[0]/100,min(ymi[0], min(mis[nonzero(mis)])/1.1))
            plt.ylim(ymi[0], yma[0])
            if hky85: plt.ylabel(r"$\theta_{f/g}$")
            else: plt.ylabel(r"$\theta_f$")

        if hky85:
            plt.tick_params('x', labelbottom=False)

        if hky85 and hky85_i:
            plt.subplot(312, sharex=ax1)
            plt.errorbar(mdepths, mthetas_r, color=col[i], xerr=sdepths, yerr=sthetas_r, fmt=fmts[i%len(fmts)], markersize=mks[i],linewidth=2, capsize=6)
            plt.tick_params('x', labelbottom=False)

            if args.relative:
                plt.yscale("linear")
                plt.ylim(0, 1.5)
                plt.ylabel(r"$\theta_r/\theta_{r0}$")
            else:
                plt.yscale("log")
                plt.ylabel(r"$\theta_r$")

                mas = mthetas_r + sthetas_r
                mis = mthetas_r - sthetas_r
                yma[1] = max(yma[1], max(mas[nonzero(mas)])*1.1)
                ymi[1] = min(ymi[1], min(mis[nonzero(mis)])/1.1)
                plt.ylim(ymi[1], yma[1])

            plt.subplot(313, sharex=ax1)
            plt.errorbar(mdepths, mmus, color=col[i], xerr=sdepths, yerr=smus, fmt=fmts[i%len(fmts)], markersize=mks[i],linewidth=2, capsize=6)

            if args.relative:
                plt.yscale("linear")
                plt.ylim(0, 1.5)
                plt.ylabel(r"$\mu/\mu_0$")
            else:
                plt.yscale("log")
                plt.ylabel(r"$\mu$")

                mas = mmus + smus
                mis = mmus - smus
                yma[2] = max(yma[2], max(mas[nonzero(mas)])*1.1)
                ymi[2] = min(ymi[2], min(mis[nonzero(mis)])/1.1, yma[2]/5)
                plt.ylim(ymi[2], yma[2])

    # All
    plt.xlabel(r"Depth")
    xxs = r_[100,50,20,10,5,2,1,0.5,0.2,0.1,0.05,0.02,0.01,0.005,0.002,0.001,0.0005,0.0002,0.0001,0.00005,0.00002,0.00001]
    if xma > 1.1:
        plt.xticks(xxs, ["%2.2f"%(x) for x in xxs])
    else:
        plt.xticks(xxs, ["%.0e"%(x) for x in xxs])

    plt.xlim(xma*1.1, xmi/1.1)

    plt.tight_layout()
    plt.show()
    #plt.savefig(args.files[0].split('_')[0] + ".png")
