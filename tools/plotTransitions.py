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
    parser = argparse.ArgumentParser(description="Plot Transition tables")
    parser.add_argument("file", help="Transitions File")
    parser.add_argument("--chr", "-c", default="All", help="chr to plot")
    parser.add_argument("--mate", "-m", default="Mate1", help="Mate to plot")
    parser.add_argument("--strand", "-s", default="Fwd", help="Strand to plot")
    parser.add_argument("--end", "-e", default="5", help="Strand to plot")
    parser.add_argument("--bases", "-b", default="ACGT", help="Strand to plot")

    args = parser.parse_args()
    prop = genfromtxt(args.file, dtype=str, usecols=(0,1,2,3), skip_header=1)

    ii = where(prop[:,0] == args.chr)[0]
    ii = ii[where(prop[ii,1] == args.mate)[0]]
    ii = ii[where(prop[ii,2] == args.strand)[0]]
    ii = ii[where(prop[ii,3] == args.end)[0]]

    data = genfromtxt(args.file, skip_header=1)[ii]

    ymax = 0.
    for i in range(5, 21):
        ymax = max(ymax, max(data[0:30, i]))


    Nbases = len(args.bases)
    cols = {"A": r_[5:9], "C" : r_[9:13], "G" : r_[13:17], "T" : r_[17:21]}

    ax = plt.subplot(111)
    for i,b in enumerate(args.bases):
        axi = plt.subplot(Nbases*100+10 + i + 1, sharex=ax)
        c   = cols[b] 
        plt.tick_params('x', labelbottom=False)
        plt.plot(data[:,4], data[:,c[0]], label="{}-A".format(b))
        plt.plot(data[:,4], data[:,c[1]], label="{}-C".format(b))
        plt.plot(data[:,4], data[:,c[2]], label="{}-G".format(b))
        plt.plot(data[:,4], data[:,c[3]], label="{}-T".format(b))
        plt.legend()
        plt.ylim([0, 1.1*ymax])

    ax.set_xlabel("position")
    ax.set_title("{}: {}, {}-Strand, {}-End".format(args.chr, args.mate, args.strand, args.end))


    plt.show()
