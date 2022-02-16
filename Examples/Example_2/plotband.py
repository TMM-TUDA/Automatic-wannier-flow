import matplotlib.pyplot as plt
import numpy as np
import re
import os
from pymatgen.io.vasp.outputs import Procar
from pymatgen.io.vasp.outputs import Outcar

plt.rc('text', usetex=False)

def readpro(filen='PROCAR'):
    data = (list(Procar(filen).data.items())[0][1])

    def energy(x): return re.match('.+energy\s+([-+]?\d*\.\d+|[-+]?\d+).+', x)

    def kbi(x): return re.match('[^\d]+(\d+)[^\d]+(\d+)[^\d]+(\d+)[^\d]+', x)

    enlist = []
    with open(filen, 'r') as f:
        (f.readline())
        kbid = kbi(f.readline())
        kp, band, ion = (int(kbid.group(1)), int(
            kbid.group(2)), int(kbid.group(3)))
        for i in f:
            m = energy(i)
            if m != None:
                enlist.append(m.group(1))
    enlist = np.float_(enlist)
    if kp * band == len(enlist):
        enlist = (enlist.reshape(kp, band))
        print(ion)
        for i in range(ion):
            with open(str(i + 1) + '.dat', 'w') as f:
                f.write('kpoint      energy         s      p      d    \n')
                for j in range(band):
                    for k in range(kp):
                        tmp = tuple(data[k][j][i])
                        f.write('%5s%12.6f%12.6f%12.6f%12.6f\n'
                                % (k + 1, enlist[k, j], tmp[0], tmp[1], tmp[2]))
                    f.write('\n')
    elif 2 * kp * band == len(enlist):
        enlistup = enlist[:kp * band]
        enlistup = (enlistup.reshape(kp, band))
        enlistdown = enlist[kp * band:]
        enlistdown = (enlistdown.reshape(kp, band))
        for i in range(ion):
            with open(str(i + 1) + '.dat.up', 'w') as f:
                f.write('kpoint      energy         s      p      d    \n')
                for j in range(band):
                    for k in range(kp):
                        tmp = tuple(data[k][j][i])
                        f.write('%5s%12.6f%12.6f%12.6f%12.6f\n'
                                % (k + 1, enlistup[k, j], tmp[0], tmp[1], tmp[2]))
                    f.write('\n')
        for i in range(ion):
            with open(str(i + 1) + '.dat.down', 'w') as f:
                f.write('kpoint      energy         s      p      d    \n')
                for j in range(band):
                    for k in range(kp):
                        tmp = tuple(data[k][j][i])
                        f.write('%5s%12.6f%12.6f%12.6f%12.6f\n'
                                % (k + 1, enlistdown[k, j], tmp[0], tmp[1], tmp[2]))
                    f.write('\n')
    else:
        raise IOError('The Procar is not correct!')
    return ion

def band(efermi=0, save='./band.png', klist='./KPOINTS', data='./1.dat',plotrange=[-5,5]):
    if not os.path.isfile(data):
        readpro('./PROCAR')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if os.path.isfile(data):
        bub = np.loadtxt('1.dat', comments='k')
    elif os.path.isfile(data + '.up') and os.path.isfile(data + '.down'):
        bub = np.loadtxt('1.dat.up', comments='k')
        bubdn = np.loadtxt('1.dat.down', comments='k')
        ydn = bubdn[:, 1]

    x = bub[:, 0]
    y = bub[:, 1]
    npoints = 0
    for i in x:
        npoints = int(i)
        if npoints > i:
            break
    nbands = len(x) // npoints
    for i in range(nbands):
        plt.plot(x[i * (npoints):(i + 1) * (npoints)],
                 y[i * (npoints):(i + 1) * (npoints)], 'b-')
        try:
            plt.plot(x[i * (npoints):(i + 1) * (npoints)],
                     ydn[i * (npoints):(i + 1) * (npoints)], 'r-')
        except:
            pass
    kfile = open(klist).readlines()

    highsymp = highsymp = lambda x: re.match('.*!\s+(.+)', x)
    highpoint = []

    for i in kfile:
        try:
            highpoint.append((highsymp(i).group(1)).replace(' ', ''))
        except:
            pass
    npointperline = npoints // (len(highpoint) // 2)
    highpoint = list(map(lambda x: r'$' + x + '$', highpoint))
    xticker = []
    for i in range(len(highpoint) // 2 + 1):
        xticker.append(i * npointperline)
        ax.axvline(i * npointperline, color='k', linestyle='-')
    ax.set_xticks(xticker)
    ax.axhline(y=efermi, color='k', linestyle='-')
    plt.ylim([efermi + plotrange[0], efermi + plotrange[1]])
    bandpoints = [highpoint[0]]
    tmp = -1
    index = 0

    def merge(x, y):
        if x == y:
            return x
        else:
            return x + '$|$' + y
    for i in range(len(highpoint) - 1):
        if i % 2 == 1:
            bandpoints.append(merge(highpoint[i], highpoint[i + 1]))
    bandpoints.append(highpoint[-1])
    ax.set_xticklabels(bandpoints)
    xm= npointperline*(len(bandpoints)-1)
    plt.xlim([0,xm])
    plt.ylabel('Energy (eV)',fontsize=15)
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.savefig(save)

efermi=Outcar('../OUTCAR').efermi
if not os.path.isfile('1.dat'):
    readpro()
band(efermi=efermi,plotrange=[-1,1])
