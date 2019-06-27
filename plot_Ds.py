#! /home/kazushi/anaconda3/bin/python
import numpy as np
import matplotlib.pyplot as plt

class DKdata:
    def load(self,fname):
        fp=open(fname,"r")
        fp.readline()
        Sr=[]
        D=[]
        K=[]
        a=[]
        beta=[];
        for row in fp:
            data=list(map(float,row.strip().split(",")));
            Sr.append(data[0])
            D.append(data[1])
            beta.append(data[2])
            K.append(data[3])
            a.append(data[4])
        self.Sr=np.array(Sr)
        self.D=np.array(D)
        self.K=np.array(K)
        self.beta=np.array(beta)
        self.a=np.array(a)

    def plot_D(self,ax,msz=6):
        line,=ax.plot(self.Sr,self.D,"-o",markersize=msz)
        return(line)
    def plot_K(self,ax,msz=6):
        line,=ax.plot(self.Sr,self.K,"-o",markersize=msz)
        return(line)
    def plot_a(self,ax,msz=6):
        line,=ax.plot(self.Sr,self.a,"-o",markersize=msz)
        return(line)
    def plot_beta(self,ax,msz=6):
        line,=ax.plot(self.Sr,self.beta,"-o",markersize=msz)
        return(line)

if __name__=="__main__":

    dk=DKdata()


    fig1=plt.figure()
    fig2=plt.figure()
    fig3=plt.figure()

    ax1=fig1.add_subplot(111)
    ax2=fig2.add_subplot(111)
    ax3=fig3.add_subplot(111)
    ax=[ax1,ax2,ax3]

    fsz=14;
    for m in range(len(ax)):
        ax[m].grid(True)
        ax[m].set_xlabel("Degree of Saturation Sr",fontsize=fsz)
        ax[m].set_xlim([0,1])
        ax[m].tick_params(labelsize=fsz)
    ax[0].set_ylabel("Diffusion constant $D$",fontsize=fsz)
    ax[1].set_ylabel("Diffusion constant $K$",fontsize=fsz)
    ax[2].set_ylabel("Exponent of time $a$",fontsize=fsz)

    nums=[5,6,3,4]
    lD=[]
    lK=[]
    la=[]
    for k in nums:
        fname="DATA"+str(k)+"/deff.dat"
        dk.load(fname)
        lD.append(dk.plot_D(ax[0],msz=8))
        lK.append(dk.plot_K(ax[1],msz=8))
        la.append(dk.plot_a(ax[2],msz=8))

    lbl=["Model 1","Model 2","Model 3","Model 4"]
    ax[0].legend(lD,lbl)
    ax[1].legend(lK,lbl)
    ax[2].legend(la,lbl)

    fig1.savefig("Ds.png",bbox_inches="tight")
    fig2.savefig("Ks.png",bbox_inches="tight")
    fig3.savefig("aexps.png",bbox_inches="tight")
    plt.show()


