import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import logit
from scipy.special import expit

def modelFn(model):
    coefs = model.split("+")
    vs = []
    es = []
    for c in coefs:
        c = c.strip()
        if c.startswith("x"):
            c = "1*" + c
        spl = c.split("*")
        vs.append(float(spl[0]))
        if len(spl) > 1:
            if spl[1] == "x":
                es.append(1)
            else:
                spl = spl[1].split("x")
                es.append(int(spl[1].split("^")[1]))
        else:
            es.append(0)

    def fn(x):
        val = 0.
        for i, v in enumerate(vs):
            val += v*(x**es[i])
        return val

    return fn
        
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot recal model(s)")
    parser.add_argument("models", nargs='*', default=["0 + x"], help="models to parse")

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
        plt.plot(qs, recals, label=args.models[i])

    plt.legend()
    plt.xlabel("quality")
    plt.ylabel("recal")
    plt.show()
