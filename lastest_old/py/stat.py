#!/usr/bin/python3
# Filename: stat.py

import numpy as np
#np.seterr(divide='ignore', invalid='ignore')
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
import setinfo

fflg=str(sys.argv[1]) #the filename to plot
clmn=int(sys.argv[2]) #the column to plot

myinfo = setinfo.Infolog('info.now')
myinfo.add('put.rc')
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

x1=int(myinfo.dk['stat-x1']);x2=int(myinfo.dk['stat-x2']);
plt.hist(y,bins=200,normed=True,range=(x1,x2), histtype='step')
y1,y2=plt.ylim()
plt.ylabel(myinfo.dk['stat-ylabel'])
plt.xlabel(statname)
plt.title(myinfo.dk['stat-title'])

ymean=y.mean()
plt.axvline(x=ymean,color='r',ls='--')
print('the mean value of '+statname+'\n', ymean)

txt = plt.text(x1+0.75*(x2-x1),0.7*y2,
'temp=%.2e\n'%myinfo.dk['temp']
+'bead=%d\n'%myinfo.dk['bead']
+'dt=%.1e\n'%myinfo.dk['dtime']
#+'gamma=%.1e\n'%myinfo.dk['gammaAD']
+'\nmean=%f'%ymean,
bbox=dict(facecolor='grey', alpha=0.05))
txt.xycoords = 'figure pixels'

plt.savefig(fflg+'_'+statname+'.png')

if(sys.argv[3]=='Y'):
	plt.show()
	#os.system('echo %.10e > .mth.tmp'%ymean)


