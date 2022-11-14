import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import logit
from scipy.special import expit

def modelFn(model):
    model = "lambda x: " + model
    fn = eval(model)
    return fn

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot recal model(s)")
    parser.add_argument("models", nargs='*', default=["x"], help="models to parse")

    args = parser.parse_args()

    fns = []
    for m in args.models:
        fns.append(modelFn(m))

    qs = np.r_[0:93:1]
    ps = 10**(-qs/10)
    vs = logit(ps)
    for i, fn in enumerate(fns):
        etas = fn(vs)
        recals = -10*np.log10(expit(etas))
        plt.plot(qs, etas, label="f" + str(i+1) + "=" + args.models[i])
        plt.plot(qs, recals, label="logit(f" + str(i+1) + ")")

    plt.legend()
    plt.xlabel("quality")
    plt.ylabel("recal")
    plt.show()
