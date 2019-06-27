#!/home/kazushi/anaconda3/bin/python

import numpy as np
import matplotlib.pyplot as plt
import gen_cmap as gcmap

class Pixel:
    def __init__(self,Wx,Wy,Nx,Ny):
        self.wx=Wx;
        self.wy=Wy;
        self.nx=Nx;
        self.ny=Ny;
        self.dx=Wx/Nx;
        self.dy=Wy/Ny;

    def draw(self,ax,ix,iy,lw=1.0):
        dx=self.dx
        dy=self.dy
        x0=ix*dx;
        y0=iy*dy;

        xs=[x0,x0+dx,x0+dx,x0,x0]
        ys=[y0,y0,y0+dy,y0+dy,y0]

        ax.plot(xs,ys,"k",linewidth=lw)

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

        self.indx=ii;
        self.jndx=jj;
        self.phs=DAT[:,1]
        self.Xa=Xa;
        self.Xb=Xb;
        self.Ndiv=Ndiv;
        self.kcell[ii,jj]=DAT[:,1]
        self.kcell=np.transpose(self.kcell)
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
        ax.tick_params(labelsize=fsz)
        if show_fig:
            plt.show()
        return(fig)
    def show_kcell_ax(self,ax):
        ext=[self.Xa[0],self.Xb[0],self.Xa[1],self.Xb[1]];
        cm = gcmap.generate_cmap(['lightyellow', 'deepskyblue', 'dimgray'])    # air, water, solid
        im=ax.imshow(self.kcell,extent=ext,interpolation="none",vmin=0,vmax=2,cmap=cm,origin="lower")
        #im=ax.imshow(self.kcell,extent=ext,interpolation="none",vmin=0,vmax=2,cmap=cm)
        fsz=14
        ax.set_xlabel("x/W",fontsize=fsz)
        ax.set_ylabel("y/W",fontsize=fsz)
        ax.set_title(self.data_file,fontsize=fsz)
        ax.tick_params(labelsize=fsz)
    def mkimg(self,fname="model.png"):
        fig=self.show_kcell(show_fig=False)
    def mkimg(self,fname="model.png"):
        fig=self.show_kcell(show_fig=False)
        fig.savefig(fname,bbox_inches="tight");
    def show_porecells(self,ax):
        px=Pixel(1.0,1.0,self.Ndiv[0],self.Ndiv[1])

        ncell=len(self.indx)
        for k in range(ncell):
            ii=self.indx[k]
            jj=self.jndx[k]
            px.draw(ax,ii,jj,lw=0.5)



if __name__=="__main__":

    pdat=PORE()

    nums=np.arange(9)+1
    nums=nums.astype(int)
    seed_num=1
    for k in nums:
        if k==9:
            fname="Sr"+str(k)+"/pore.dat"
        else:
            fname="Sr"+str(k)+"/SEED"+str(seed_num)+"/pore.dat"
        pdat.load(fname)
        fimg="sr"+str(k)+".png"
        pdat.mkimg(fname=fimg)

    kcell=pdat.kcell/2
    pdat.kcell=kcell.astype(int)*2
    pdat.mkimg("sr0.png")




