import re,os,sys
from time import time
import numpy as np
from numpy.linalg import eigvals
from numpy.linalg import inv
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.patches as mpatches
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
from matplotlib.colors import LinearSegmentedColormap
from multiprocessing import Pool
from pymatgen.io.vasp.outputs import Outcar

## Author: Zhang Zeying

plt.rc('text', usetex=False)

def wband(efermi=0, bandf='band/1.dat', klist='band/KPOINTS',plotrange=[8,5]):
    '''
    this part is for comparing the wannier functions
    '''
    fig = plt.figure()
    ax = fig.add_subplot(111)
    bub = np.loadtxt(bandf, comments='k')
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
    plt.ylim([efermi - plotrange[0], efermi + plotrange[1]])
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

class readhr(object):
    def __init__(self):
        filename = 'wann/wannier90_hr.dat'
        hrdat = open(filename, 'r')
        hrdat.readline()
        self.nwann = int(hrdat.readline())
        self.nrpts = int(hrdat.readline())
        self.dengen = []
        if self.nrpts % 15 == 0:
            ndeline = self.nrpts // 15
        else:
            ndeline = self.nrpts // 15 + 1
        for i in range(ndeline):
            deg = hrdat.readline()
            deg = str.split(deg)
            for j in deg:
                self.dengen.append(int(j))
        self.hmnr = np.zeros(
            (self.nrpts, self.nwann, self.nwann), dtype=complex)
        self.ws = np.zeros((self.nrpts, 3), dtype=int)
        for i in range(self.nrpts):
            for j in range(self.nwann):
                for k in range(self.nwann):
                    hsigleline = str.split(hrdat.readline())
                    tem1 = float(hsigleline[-2])
                    tem2 = float(hsigleline[-1])
                  #  self.hmnr[i, k, j] = (tem1 + 1j * tem2)
                    self.hmnr[i, k, j] = (tem1 + 1j * tem2) / self.dengen[i]
            self.ws[i, 0] = int(hsigleline[0])
            self.ws[i, 1] = int(hsigleline[1])
            self.ws[i, 2] = int(hsigleline[2])
            self.ws = np.array(self.ws)
        hrdat.close()

    def hk(self, k):
        hk = np.zeros((self.nwann, self.nwann))
        for i in range(self.nrpts):
            hk = hk + (self.hmnr[i, :, :] * np.exp(2.0 *
                                                   np.pi * 1j * np.inner(k, self.ws[i])))
        return hk

    def readkp(self, filen='band/KPOINTS'):
        def highsymp(x): return re.match('.*!\s+(.+)', x)

        def highsymk(x): return re.match(
            '([-+]?\d*\.\d+\s+|[-+]?\d+\s+){3}', x)
        highpoint = []
        highk = []
        kfile = open(filen).readlines()
        pointsperline = int(kfile[1])
        for i in kfile:
            symp = highsymp(i)
            symk = highsymk(i)
            if symk != None:
                highk.append(np.float_(str.split(symk.group())))
            if symp != None:
                highpoint.append((highsymp(i).group(1)).replace(' ', ''))

        def kline(start, end, point):
            x = np.linspace(start[0], end[0], point)
            y = np.linspace(start[1], end[1], point)
            z = np.linspace(start[2], end[2], point)
            return np.transpose(np.array([x, y, z]))

        klist = []
        for i in range(len(highk) // 2):
            klist.append(kline(highk[2 * i], highk[2 * i + 1], pointsperline))
        klist = np.array(klist).reshape((pointsperline) * len(highk) // 2, 3)
        return klist
    def nocc(self, efermi=0.0):
        klist=[[0,0,0],[0,0,0.5],[0,0.5,0],[0.5,0,0],\
                [0.5,0.5,0],[0.5,0,0.5],[0,0.5,0.5],[0.5,0.5,0.5],[0.2,0.3,0.4]]
        band=[]
        bandt=[]
        for i in klist:
            band.append(np.sort(eigvals(self.hk(i)).real)-efermi)
            bandt.append(np.sort(eigvals(self.hk(i)).real))
        occs=[]
        for i in band:
            num=0
            for j in i:
                if j<0:
                    num=num+1
                else:
                    occs.append(num)
                    break
        mean=np.mean(occs)
        dec=int(mean)-mean
        if dec != 0:
            print('WARNING!! OCC NUMBER MAY NOT CORRECT! OCC FOR EACH POINT IS',occs)
        return int(np.round(mean))

    def bandplot(self, efermi=0.0):
        band = []
        for i in self.readkp():
            band.append(np.sort(eigvals(self.hk(i)).real))
        band = (np.transpose(band))
        x = np.arange(len(band[1])) + 1
        for i in band:
            plt.plot(x, i, '-r')
        plt.show()

    def combandplot(self, efermi=0.0, save='./', bandf='band/1.dat', klist='band/KPOINTS',winf='wann/wannier90.win', plotrange=[8,5]):
        band = []
        for i in self.readkp():
            band.append(np.sort(eigvals(self.hk(i)).real))
        band = (np.transpose(band))
        x = np.arange(len(band[1])) + 1
        wband(efermi=efermi, bandf=bandf, klist=klist, plotrange=plotrange)
        for i in band:
            plt.plot(x, i, '--r')
        red_patch = mpatches.Patch(color='red', label='Wannier')
        blue_patch = mpatches.Patch(color='blue', label='DFT')
        plt.legend(handles=[red_patch, blue_patch],loc='upper right')
        LSQ=self.comparedata(bandf=bandf,winf=winf)
        plt.title('Least Squares=%12.6f'%(LSQ))
        plt.ylabel('Energy (eV)',fontsize=15)
        plt.tick_params(axis='both', which='major', labelsize=15)
        plt.savefig(save + 'wannband')
        plt.close()
        return LSQ

    def comfw(self, efermi=0.0,eig='../band/EIGENVAL',klist='../band/KPOINTS',save='./wann',pr=[-1,1],LSQ=False,winf='./wannier90.win',eig1='../band/EIGENVAL',eig2='../dos/EIGENVAL'):
        bandp=band.bandstr()
        bandp.band(filen=eig)
        bandp.readk(kp=klist)
        xticker = []
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.axhline(y=efermi, color='k', linestyle='-')
        plt.ylim([efermi +pr[0], efermi + pr[1]])
        xm= bandp.npointperline*(len(bandp.bandpoints)-1)
        plt.xlim([0,xm])
        for i in range(len(bandp.bandpoints)):
            xticker.append(i * bandp.npointperline)
            ax.axvline(i * bandp.npointperline, color='k', linestyle='-')
        ax.set_xticklabels(bandp.bandpoints)
        ax.set_xticks(xticker)
        if (bandp.spin)=='1':
            for i in range(bandp.nband):
                plt.plot(bandp.bdup[:,i],'b-')
        elif bandp.spin=='2':
            for i in range(bandp.nband):
                plt.plot(bandp.bdup[:,i],'b-')
                plt.plot(bandp.bddn[:,i],'r-')
        bandfs = []
        for i in self.readkp(filen=klist):
            bandfs.append(np.sort(eigvals(self.hk(i)).real))
        bandfs = (np.transpose(bandfs))
        x = np.arange(len(bandfs[1])) 
        for i in bandfs:
            plt.plot(x, i, '--r')
        red_patch = mpatches.Patch(color='red', label='Wannier')
        blue_patch = mpatches.Patch(color='blue', label='DFT')
        plt.legend(handles=[red_patch, blue_patch],loc='upper right')
        if LSQ:
            lsq1=self.comeig(eig=eig1,winf=winf,efermi=efermi)
            if eig2 is not None:
                lsq2=self.comeig(eig=eig2,winf=winf,efermi=efermi)
                plt.title('Least Squares=Band%12.6f Sactter%12.6f'%(lsq1,lsq2))
            else:
                plt.title('Least Squares=Band%12.6f'%(lsq1))
        else:
            pass
        plt.ylabel('Energy (eV)',fontsize=15)
        plt.tick_params(axis='both', which='major', labelsize=15)
        if isinstance(save, str):
            plt.savefig(save)
        else:
            plt.show()
        plt.close()

    def printdata(self):
        klist = self.readkp('../KPOINTS')
        band = []
        for i in klist:
            band.append(np.sort(eigvals(self.hk(i)).real))
        band = (np.transpose(band))
        eig = open('eig', 'w')
        for i in band:
            for j in i:
                eig.write('%.8E\n' % j)
            eig.write('\n')

    def comeig(self,eig='../band/EIGENVAL',winf='./wannier90.win',efermi=0):
        def diswin(x):
            return re.match('dis_win_min=([-+]?\d*\.\d+\s+|[-+]?\d+\s+)', x)
        def lsq(list1,list2):
            s=0
            for i,j in zip(list1,list2):
                s=s+(i-j)*(i-j)
            return s
        with open(winf,'r') as f:
            for i in f:
                if (diswin(i)) is not None:
                    dw=float(diswin(i).group(1))
                    break
        bandp=band.bandstr()
        bandp.band(filen=eig)
        start=0
        for i in range(bandp.nband):
            if min(bandp.bdup[:,i])<dw:
                start+=1
            else:
                break
        wannband=[]
        for i in range(bandp.nk):
            wannband.append(np.sort(eigvals(self.hk(bandp.klist[:,i])).real))
        wannband = (np.transpose(wannband))
        sums=0
        cutband=0
        for i in range(self.nwann):
            if max(wannband[i]) < efermi:
                cutband+=1
            else:
                break
        for i in range(cutband):
            sums=sums+lsq(wannband[i],bandp.bdup[:,start+i])
        npoints=len(bandp.bdup[:,0])
        res=100*sums/(bandp.nband*npoints)
        return res

    def comparedata(self,winf='wann/wannier90.win',bandf='band/1.dat'):
        def diswin(x): 
            return re.match('dis_win_min=([-+]?\d*\.\d+\s+|[-+]?\d+\s+)', x)
        def lsq(list1,list2):
            s=0
            for i,j in zip(list1,list2):
                s=s+(i-j)*(i-j)
            return s
        with open(winf,'r') as f:
            for i in f:
                if (diswin(i)) is not None:
                    dw=float(diswin(i).group(1))
                    break
        bub = np.loadtxt(bandf, comments='k')
        x = bub[:, 0]
        y = bub[:, 1]
        npoints = 0
        for i in x:
            npoints = int(i)
            if npoints > i:
                break
        nbands = len(x) // npoints
        bands=np.zeros((nbands,npoints))
        for i in  range(nbands):
            bands[i,:]=y[i * (npoints):(i + 1) * (npoints)]
        start=0
        for i in bands:
            if min(i)<dw:
                start+=1
            else:
                break
        wannband = []
        for i in self.readkp('band/KPOINTS'):
            wannband.append(np.sort(eigvals(self.hk(i)).real))
        wannband = (np.transpose(wannband))
        sums=0
        for i in range(self.nwann//2):
            sums=sums+lsq(wannband[i],bands[start+i])
        return 100*sums/(nbands*npoints)

if os.path.isfile('wann/wannier90_hr.dat') and os.path.isfile('band/PROCAR') :
    efermi=Outcar('OUTCAR').efermi
    hr=readhr()
    readhr().combandplot(efermi=efermi,plotrange=[2,2])
else:
    print("Files not found")
