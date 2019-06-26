#!/home/kazushi/anaconda3/bin/python
import numpy as np
import matplotlib.pyplot as plt

class Sct:
    def load(self,fname):
        fp=open(fname,"r")
        x=[]; y=[];
        ux=[]; uy=[];
        fp.readline()
        for var in fp:
            dat=var.strip().split(" ")
            x.append(float(dat[0]));
            y.append(float(dat[1]));
            ux.append(float(dat[2]));
            uy.append(float(dat[3]));

        self.x=np.array(x);
        self.y=np.array(y);
        self.ux=np.array(ux);
        self.uy=np.array(uy);
    def plot_xy(self,ax,msz=1):
        ax.plot(self.x,self.y,".g",markersize=msz)
    def plot_uv(self,ax,msz=1,clr="g"):
        ax.plot(self.ux,self.uy,"."+clr,markersize=msz)

    def uhist(self,ax,nbin=30,clr="b"):
        u=np.hstack([self.ux,self.uy])
        ax.hist(u,bins=nbin,normed=True,color=clr)


if __name__=="__main__":


    fig1=plt.figure()
    fig2=plt.figure()

    ax1=fig1.add_subplot(111)
    ax2=fig2.add_subplot(111)
    ax=[ax1,ax2]


    bins=50
    pts=Sct()
    nums=[1,2,3,4]
    clrs=["b","g","m","r","y"]
    nclrs=len(clrs)
    H_max=3


    for k in nums:
        ax[0].set_ylim([-3,3])
        ax[0].set_aspect(1.0)
        ax[0].grid(True)
        ax[0].set_xlim([-3,3])

        ax[1].set_ylim([0,H_max])
        ax[1].set_xlim([-3,3])
        ax[1].grid(True)

        fname="rwk"+str(k)+".out"
        pts.load(fname)
        pts.plot_uv(ax[0],clr="g")
        hst=pts.uhist(ax2,nbin=bins,clr="g")

        fsct="sctt"+str(k)+".png"
        fhist="hist"+str(k)+".png"
        ax1.set_title(fname)
        ax2.set_title(fname)
        ax1.tick_params(labelsize=14)
        ax2.tick_params(labelsize=14)

        ax1.set_xlabel("x/W",fontsize=14)
        ax1.set_ylabel("y/W",fontsize=14)
        ax2.set_xlabel("u$_i$/W",fontsize=14)
        ax2.set_ylabel("Probability density",fontsize=14)

        fig1.savefig(fsct,bbox="inches_tight")
        fig2.savefig(fhist,bbox="inches_tight")

        ax1.clear()
        ax2.clear()


