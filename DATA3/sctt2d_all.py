#!/home/kazushi/anaconda3/bin/python
import numpy as np
import matplotlib.pyplot as plt

class Sct:
    def __init__(self):
        self.x=np.array([])
        self.y=np.array([])
        self.ux=np.array([])
        self.uy=np.array([])
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
        ax.hist(u,bins=nbin,normed=True,color=clr,log=True)


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


    PTS=Sct()
    for k in nums:  # RWK step
        ax[0].set_ylim([-3,3])
        ax[0].set_aspect(1.0)
        ax[0].grid(True)
        ax[0].set_xlim([-3,3])

        #ax[1].set_ylim([0,1.2])
        ax[1].set_xlim([-3,3])
        ax[1].grid(True)

        Dir="SEED"
        fname="rwk"+str(k)+".out"

        for m in range(10):
            Fname=Dir+str(m+1)+"/"+fname
            pts.load(Fname)
            PTS.x=np.hstack([PTS.x,pts.x])
            PTS.y=np.hstack([PTS.y,pts.y])
            PTS.ux=np.hstack([PTS.ux,pts.ux])
            PTS.uy=np.hstack([PTS.uy,pts.uy])

        sig=np.std(np.hstack([PTS.ux,PTS.uy]))
        uu=np.linspace(-3,3,201)
        G=np.exp(-uu*uu/(2.*sig*sig))/sig/np.sqrt(2*np.pi)

        PTS.plot_uv(ax[0],clr="g")
        hst=PTS.uhist(ax2,nbin=bins,clr="g")
        ax2.semilogy(uu,G,"b",linewidth=2)
        ax2.set_ylim([1.e-04,10])

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


