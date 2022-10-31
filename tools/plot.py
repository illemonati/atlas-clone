import argparse
from numpy import *
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

    xs = arange(0, 100, 0.01)
    for i, fn in enumerate(fns):
        plt.plot(xs, fn(xs), label=args.models[i])

    plt.legend()
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.show()
