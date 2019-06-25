#! /home/kazushi/anaconda3/bin/python
import numpy as np
import matplotlib.pyplot as plt
import gen_cmap as gcmap


fname="kcell3.dat";

nums=np.arange(19);


#cm = gcmap.generate_cmap(['lightgray', 'deepskyblue', 'saddlebrown'])    # air, water, solid
cm = gcmap.generate_cmap(['lightyellow', 'deepskyblue', 'dimgray'])    # air, water, solid
for k in nums:
    fname="kcell"+str(k)+".dat"
    fn_out="model"+str(k)+".png"
    fp=open(fname,"r");

    fp.readline(); # computational domain
    Xa=list(map(float,fp.readline().lstrip().split(",")));
    Xb=list(map(float,fp.readline().lstrip().split(",")));

    fp.readline(); # computational domain
    Ndiv=list(map(int,fp.readline().lstrip().split(",")));

    fp.readline();	# Imaging area
    K=list(map(float,fp.readlines()));	# Imaging area
    K=np.transpose(np.reshape(K,Ndiv))

    print(np.shape(K));
    fig=plt.figure();
    ax=fig.add_subplot(111)
    im=ax.imshow(K,origin="lower",extent=[Xa[0],Xb[0],Xa[1],Xb[1]],interpolation="none",vmin=0,vmax=2,cmap=cm);
    #plt.colorbar(im);
    #plt.show();
    fig.savefig(fn_out,bbox_inches="tight")


fp.close();
