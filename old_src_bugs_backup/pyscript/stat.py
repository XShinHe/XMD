#!/usr/bin/python3
# Filename: stat.py

import numpy as np
#np.seterr(divide='ignore', invalid='ignore')
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
import Cata

fflg=str(sys.argv[1]) #the filename to plot
clmn=int(sys.argv[2]) #the column to plot

cata = Cata.Catalog('_plot.in')
cata.add('_put.in')
#

a=pd.read_csv(fflg,sep='\s+')
statname=a.columns.values[clmn]
mtx=a.values.T
if(clmn>0):
    y=mtx[clmn:clmn+1,]
elif(clmn==0):
    y_all=mtx[1:,]
    t=y_all.reshape(1,-1)
else:
    pass
y=y[0]

x1=int(cata.dk['x1']);x2=int(cata.dk['x2']);
plt.hist(y,bins=200,normed=True,range=(x1,x2), histtype='step')
y1,y2=plt.ylim()
plt.ylabel(cata.dk['ylabel'])
plt.xlabel(cata.dk['xlabel'])
plt.title(cata.dk['title'])

ymean=y.mean()
plt.axvline(x=ymean,color='r',ls='--')
print('the mean value of '+statname+'\n', ymean)

plt.text(0.7*x2,0.7*y2,'bead=%d\n'%cata.dk['bead']+'dt=%f\n'%cata.dk['dtime']+'temp=%f\n'%cata.dk['temp']+'\n\n\nmean=%f'%ymean)
plt.savefig(fflg+'%d'%clmn+'.png')
#plt.show()
os.system('echo %.10e > .mth.tmp'%ymean)




