#!/home/kazushi/anaconda3/bin/python
import numpy as np

import matplotlib.pyplot as plt

fname="Sr9/rwk.out"
fp=open(fname,"r")

x=[];
y=[];
ux=[];
uy=[];
fp.readline()
for var in fp:
    dat=var.strip().split(" ")
    x.append(float(dat[0]));
    y.append(float(dat[1]));
    ux.append(float(dat[2]));
    uy.append(float(dat[3]));


fig=plt.figure()
ax1=fig.add_subplot(211)
ax2=fig.add_subplot(212)
#ax.plot(x,y,".")
ax1.hist(ux,bins=30);
ax2.hist(uy,bins=30)
ax1.grid(True)
ax2.grid(True)

plt.show()



