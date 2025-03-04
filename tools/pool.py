from numpy import *

if __name__ == "__main__":
    Ns = [0, 0, 10, 100, 40, 100, 10, 10, 10, 10, 100, 40, 0, 30, 10, 0, 100, 30, 30, 30, 100, 40, 30, 30]
    Nmin = 50
    print(Ns)

    iis = []
    nns = []
    for i, n in enumerate(Ns):
        if n > 0:
            nns.append(n)
            iis.append([i])

    print(iis)
    print(nns)

    # lower
    while nns[0] < Nmin:
        nns[1] += nns[0]
        iis[1].extend(iis[0])
        del nns[0]
        del iis[0]

    # Upper
    while nns[-1] < Nmin:
        nns[-2] += nns[-1]
        iis[-2].extend(iis[-1])
        del nns[-1]
        del iis[-1]
        

    iMin = argmin(nns)
    while nns[iMin] < Nmin:
        if nns[iMin - 1] < nns[iMin + 1]:
            nns[iMin - 1] += nns[iMin]
            iis[iMin - 1].extend(iis[iMin])
        else:
            nns[iMin + 1] += nns[iMin]
            iis[iMin + 1].extend(iis[iMin])
        del nns[iMin]
        del iis[iMin]
        iMin = argmin(nns)

    for i, ii in enumerate(iis):
        print(ii, i, nns[i])
