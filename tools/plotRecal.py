import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import logit
from scipy.special import expit

def T(x):
    return logit(10**(-x/10))

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

    xs = np.r_[0:150:1]
    for i, fn in enumerate(fns):
        eta = fn(xs)
        ax1 = plt.subplot(311)
        plt.plot(xs, eta, label="f" + str(i+1) + " = " + args.models[i])
        plt.tick_params('x', labelbottom=False)
        plt.ylabel("eta")
        plt.legend()

        prob = expit(eta)
        plt.subplot(312, sharex=ax1)
        plt.tick_params('x', labelbottom=False)
        plt.plot(xs, prob)
        plt.ylabel("prob")

        plt.subplot(313, sharex=ax1)
        q = -10*np.log10(prob)
        plt.plot(xs, q)
        plt.ylabel("new quality")

    plt.xlabel("x")
    plt.show()
