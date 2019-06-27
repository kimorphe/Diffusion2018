#! /home/kazushi/anaconda3/bin/python

import numpy as np
import matplotlib.pyplot as plt
# colormapをカスタマイズする
from matplotlib.colors import LinearSegmentedColormap

def generate_cmap(colors):
    """自分で定義したカラーマップを返す"""
    values = range(len(colors))
    vmax = np.ceil(np.max(values))
    color_list = []
    for v, c in zip(values, colors):
        color_list.append( ( v/ vmax, c) )
    return LinearSegmentedColormap.from_list('custom_cmap', color_list)


if __name__=="__main__":
    # 色をカスタマイズ, その1:色名称で指定
    # 参考： http://matplotlib.org/examples/color/named_colors.html
    x=np.linspace(0,2*np.pi,100)
    y=np.sin(x);
    #unique_value = set(iris.target)
    #print unique_value
    # --> [0, 1, 2]

    print(x)
    print(y)

    cm = generate_cmap(['mediumblue', 'limegreen', 'orangered'])

    fig = plt.figure(figsize=(13,7))
    im = plt.scatter(x,y, c=y, linewidths=0, alpha=.8, cmap=cm)
    fig.colorbar(im)
    plt.show()
