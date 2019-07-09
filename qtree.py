import numpy as np
import matplotlib.pyplot as plt

fname="qtree_bnd.dat"

fp=open(fname,"r");

x=[];
y=[];
for row in fp:
    dat=row.strip().split(",");
    x.append(float(dat[0]));
    y.append(float(dat[1]));

x=np.array(x);
y=np.array(y);

npx=int(len(x)/5);
x=np.reshape(x,[npx,5]);
y=np.reshape(y,[npx,5]);

x=np.transpose(x);
y=np.transpose(y);

fig,ax=plt.subplots(1,1)
ax.set_aspect(1.0);
ax.set_xlim([0.6,0.75])
ax.set_ylim([0.6,0.75])
ax.grid(True);
ax.set_xlabel("x/W",fontsize=12);
ax.set_ylabel("y/W",fontsize=12);
ax.tick_params(labelsize=12);
ax.plot(x,y,"b-",linewidth=0.5);
fig.savefig("bgrid.png",box_inches="tight");
plt.show()
