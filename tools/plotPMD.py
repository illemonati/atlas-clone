import argparse
from numpy import *
import matplotlib
import matplotlib.pyplot as plt
from scipy.special import logit
from scipy.special import expit
matplotlib.rcParams.update({'font.size': 10})

def modelFn(model):
    model = "lambda x: " + model
    fn = eval(model)
    return fn

def parseEmp(emp):
    vs = []
    spl = emp.split(",")
    for s in spl:
        vs.append(float(s))
    return r_[vs]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot recal model(s)")
    parser.add_argument("models", nargs='*', help="models to parse")
    parser.add_argument("--empiric", "-e", nargs='*', help="empiric values to parse")

    args = parser.parse_args()

    fns = []
    for m in args.models:
        fns.append(modelFn(m))

    xs = arange(0, 100, 0.01)
    for i, fn in enumerate(fns):
        plt.plot(xs, fn(xs), label=args.models[i])

    es = []
    for e in args.empiric:
        es.append(parseEmp(e))

    for i, e in enumerate(es):
        xs = arange(0, len(e))
        plt.plot(xs, e, label=str(i))


    plt.legend()
    plt.xlabel("Distance from 5'-end of Read")
    plt.ylabel("PMD C-T")
    plt.show()
