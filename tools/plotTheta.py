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
    parser.add_argument("--het",  action="store_true")

    args     = parser.parse_args()
    nFiles   = len(args.files)
    nSamples = int(ceil(nFiles/args.num))

    if (nSamples) > len(col): sys.exit("Can only compare a maximum of %d samples!"%(len(col)))

    ymi = [10., 10., 10., 10.]
    yma = [0., 0., 0., 0.]
    xmi = 1000.
    xma = 0.

    hky85 = False #at least one file is hky85
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
        ihets     = []
        ithetas_g = []
        ithetas_r = []
        ithetas_f = []
        iPMD5s    = []
        iPMD3s    = []

        for j, key in enumerate(f.readline().split()):
            key = key.decode()
            if "depth" in key: idepths.append(j)
            if "mu" in key: imus.append(j)
            if key.endswith("het"): ihets.append(j)
            if "theta_g" in key: ithetas_g.append(j)
            if "theta_r" in key: ithetas_r.append(j)
            if "thetaMLE" in key: ithetas_f.append(j)
            if "expHetMLE" in key: ihets.append(j)
            if "PMD5" in key: iPMD5s.append(j)
            if "PMD3" in key: iPMD3s.append(j)
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

        hky85_i = len(ithetas_f) == 0

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

            PMD5s = []
            for j in iPMD5s:
                PMD5s.append(data[:, j])

            PMD3s = []
            for j in iPMD3s:
                PMD3s.append(data[:, j])
        else:
            label += "_Fels"
            if args.relative:
                theta0 = data[0,ithetas_f[0]]
            else:
                theta0 = 1.

            thetas_g = []
            for j in ithetas_f:
                thetas_g.append(data[:, j]/theta0)

        fmts = ["o", "s", "X", "d", "p", "<", "^", ">"]
        lins = ["-", ":", "--"]
        mks  = [i for i in range(nFiles, 0, -1)]

        if args.median:
            print("using median")
            mdepths   = r_[[nanmedian(d) for d in depths]]
            mthetas_g = r_[[nanmedian(t) for t in thetas_g]]
            mhets     = r_[[nanmedian(h) for h in hets]]
            if hky85_i:
                mthetas_r = r_[[nanmedian(t) for t in thetas_r]]
                mmus      = r_[[nanmedian(m) for m in mus]]
                mPMD5s    = r_[[nanmedian(p) for p in PMD5s]]
                mPMD3s    = r_[[nanmedian(p) for p in PMD3s]]
        else:
            print("using mean")
            mdepths   = r_[[nanmean(d) for d in depths]]
            mthetas_g = r_[[nanmean(t) for t in thetas_g]]
            mhets     = r_[[nanmean(h) for h in hets]]
            if hky85_i:
                mthetas_r = r_[[nanmean(t) for t in thetas_r]]
                mmus      = r_[[nanmean(m) for m in mus]]
                mPMD5s    = r_[[nanmean(p) for p in PMD5s]]
                mPMD3s    = r_[[nanmean(p) for p in PMD3s]]

        sdepths   = r_[[nanstd(d) for d in depths]]
        sthetas_g = r_[[nanstd(t) for t in thetas_g]]
        shets     = r_[[nanstd(h) for h in hets]]
        if hky85_i:
            sthetas_r = r_[[nanstd(t) for t in thetas_r]]
            smus      = r_[[nanstd(m) for m in mus]]
            sPMD5s    = r_[[nanstd(p) for p in PMD5s]]
            sPMD3s    = r_[[nanstd(p) for p in PMD3s]]
        
        if hky85_i:
            print(label)
            print("depth:", mdepths)
            print("theta_g:", mthetas_g)
            print("theta_r:", mthetas_r)
            print("mu:", mmus)
            if len(mPMD5s) > 0:
                print("PMD5:", mPMD5s)
            if len(mPMD3s) > 0:
                print("PMD3:", mPMD3s)
        else:
            print(label)
            print("depth:", mdepths)
            print("theta_f:", mthetas_g)

        xmi = min(min(mdepths), xmi)
        xma = max(max(mdepths), xma)

        if hky85:
            ax1 = plt.subplot(311)
            plt.tick_params('x', labelbottom=False)
        else:
            ax1 = plt.subplot(111)

        plt.errorbar(mdepths, mthetas_g, color=col[i%nSamples], yerr=sthetas_g, fmt=fmts[i%nSamples] + lins[int(i/nSamples)], markersize=mks[i],linewidth=2, capsize=6, label=label)
        plt.hlines(mthetas_g[0], 0, 1.5*max(mdepths), col[i%nSamples], "dashed")
        plt.xscale("log")
        plt.legend(ncols=2, borderaxespad=0.)
        if args.relative:
            plt.yscale("linear")
            yma[0] = max(yma[0], max(mthetas_g[nonzero(mthetas_g)]))
            plt.ylim(0, 1.1*yma[0])
            if hky85: plt.ylabel(r"$\theta_{f/g}/\theta_0$")
            else: plt.ylabel(r"$\theta_f/\theta_0$")

        else:
            plt.yscale("log")
            yma[0] = max(yma[0], max(mthetas_g[nonzero(mthetas_g)]))
            ymi[0] = min(ymi[0], min(mthetas_g[nonzero(mthetas_g)]))
            plt.ylim(min(ymi[0]/2, yma[0]/10), yma[0]*2)
            if hky85: plt.ylabel(r"$\theta_{g}$")
            else: plt.ylabel(r"$\theta_f$")

        if hky85 and hky85_i:
            plt.subplot(312, sharex=ax1)
            plt.tick_params('x', labelbottom=False)
            plt.errorbar(mdepths, mthetas_r, color=col[i%nSamples], yerr=sthetas_r, fmt=fmts[i%nSamples] + lins[int(i/nSamples)], markersize=mks[i],linewidth=2, capsize=6)

            if args.relative:
                plt.yscale("linear")
                plt.ylim(0, 1.5)
                plt.ylabel(r"$\theta_r/\theta_{r0}$")

            else:
                plt.yscale("log")
                plt.ylabel(r"$\theta_r$")

                yma[1] = max(yma[1], max(mthetas_r[nonzero(mthetas_r)]))
                ymi[1] = min(ymi[1], min(mthetas_r[nonzero(mthetas_r)]))
                plt.ylim(min(ymi[1]/2, yma[1]/10), yma[1]*2)

        if hky85 and hky85_i:
            plt.subplot(313, sharex=ax1)

            plt.errorbar(mdepths, mmus, color=col[i%nSamples], yerr=smus, fmt=fmts[i%nSamples] + lins[int(i/nSamples)], markersize=mks[i],linewidth=2, capsize=6)
            plt.yscale("log")
            yma[2] = max(yma[2], max(mmus[nonzero(mmus)]))
            ymi[2] = min(ymi[2], min(mmus[nonzero(mmus)]))
            plt.ylim(ymi[2]/2, yma[2]*2)

            plt.ylabel(r"$\mu$")

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
