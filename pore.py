#!/home/kazushi/anaconda3/bin/python

import numpy as np
import matplotlib.pyplot as plt
import gen_cmap as gcmap


class PORE:
    def load(self,fname):
        fp=open(fname,"r")

        fp.readline()
        thE=float(fp.readline())

        fp.readline()
        Etot=float(fp.readline())

        fp.readline()
        Sr=float(fp.readline())

        fp.readline()
        tmp=fp.readline();
        Xa=list(map(float,tmp.strip().split(",")))

        tmp=fp.readline();
        Xb=list(map(float,tmp.strip().split(",")))

        fp.readline()
        tmp=fp.readline();
        Ndiv=list(map(int,tmp.strip().split(",")))

        fp.readline()
        ncell=int(fp.readline().strip())

        fp.readline()

        DAT=fp.read().strip()
        DAT=DAT.replace("\n"," ");
        DAT=DAT.split(" ")
        DAT=list(map(int,DAT))
        ndat=len(DAT)
        DAT=np.array(DAT);
        DAT=np.reshape(DAT,[int(ndat/3),3])

        self.kcell=np.ones(Ndiv)*2

        ii=DAT[:,0]/Ndiv[1];
        ii=ii.astype(int)
        jj=DAT[:,0]%Ndiv[1];

        self.Xa=Xa;
        self.Xb=Xb;
        self.Ndiv=Ndiv;
        self.kcell[ii,jj]=DAT[:,1]
        self.data_file=fname

        fp.close()
    def show_kcell(self,show_fig=True):
        fig=plt.figure()
        ax=fig.add_subplot(111)
        ext=[self.Xa[0],self.Xb[0],self.Xa[1],self.Xb[1]];
        cm = gcmap.generate_cmap(['lightyellow', 'deepskyblue', 'dimgray'])    # air, water, solid
        im=ax.imshow(self.kcell,extent=ext,interpolation="none",vmin=0,vmax=2,cmap=cm,origin="lower")
        fsz=14
        ax.set_xlabel("x/W",fontsize=fsz)
        ax.set_ylabel("y/W",fontsize=fsz)
        ax.set_title(self.data_file,fontsize=fsz)
        if show_fig:
            plt.show()
        return(fig)
    def mkimg(self,fname="model.png"):
        fig=self.show_kcell(show_fig=False)
        fig.savefig(fname,bbox_inches="tight");

if __name__=="__main__":

    pdat=PORE()

    Sr_num=5
    nums=np.arange(9)+1

    nums=nums.astype(int)
    for k in nums:
        fname="Sr"+str(Sr_num)+"/SEED"+str(k)+"/pore.dat"
        pdat.load(fname)
        #pdat.show_kcell()
        fimg="model"+str(k)+".png"
        pdat.mkimg(fname=fimg)

