# -*- coding: utf-8 -*-
"""
Created on Fri Aug 03 18:28:43 2018
Edited on Tue Oct 15 12:37:14 2019

@author: agharag
"""


from tkinter.filedialog import askopenfilename
from netCDF4 import Dataset
import numpy as np
import time
from scipy import spatial
from operator import itemgetter


''' read mesh and element conn '''

print("Select The Mesh File (.nc)")
mesh = askopenfilename() # show an "Open" dialog box and return the path to the selected file
var=input('Enter the variable name to pull: (zeta,hs,tps,dir,etc) ')

t0=time.time()


ncf=Dataset(mesh,'r')
elem  = ncf.variables['element'][:]
xnode = ncf.variables['x'][:]
ynode = ncf.variables['y'][:]
var = ncf.variables[var][:]
N_ele = elem.shape[0]
nt=var.shape[0]
N_nod= len(xnode)
nodes=np.column_stack((xnode,ynode))


''' read 63 file (maxele, water level, Hs, Tp, Dir, WindUV) '''

## read station list
print("Select The Station File (long lat)")
station = askopenfilename()
stat_list=[]
with open(station,'r') as f_st:
    f_st.readline()
    for line in f_st:
        m=line.split()
        stat_list.append([float(m[0]),float(m[1])])


''' find nearest triangles'''
t2=time.time()

### nearest k nodes to each station in stat_list 
CKD=list(spatial.cKDTree(nodes).query(stat_list,k=4)[1])


triangles=[]
for i in range(len(stat_list)): #for each station
    tri=[]    
    ps=CKD[i]                #nearest k nodes 
    for node in ps:
        tri.extend(list(np.where(elem == node+1)[0]))  # find nearest elements/triangles
    triangles.append(set(tri))    

t3=time.time()
print(' ***** found nearby elements in ', round(t3-t2,3),'s')    


## find exact triangle- and then interpolate and get data

GAMMA=[]
ADJ=[]
for i in range(len(stat_list)):
    xp,yp=stat_list[i]
    for j in triangles[i]:
        adjc=elem[j]-1
        XX  = [xp,xnode[adjc[0]],xnode[adjc[1]],xnode[adjc[2]]]
        YY  = [yp,ynode[adjc[0]],ynode[adjc[1]],ynode[adjc[2]]]

        def Area(X,Y):
            return round(abs((X[1]*Y[2]-X[2]*Y[1])-(X[0]*Y[2]-X[2]*Y[0])+(X[0]*Y[1]-X[1]*Y[0])),8)
        def Sub_Area(a,b,c):
            return Area(list(itemgetter(a,b,c)(XX)),list(itemgetter(a,b,c)(YY)))

        Tot_A  = Area(XX[1:],YY[1:])
        G0 = Sub_Area(0,2,3)/Tot_A
        G1 = Sub_Area(1,0,3)/Tot_A
        G2 = Sub_Area(1,2,0)/Tot_A
        
        if round(G0+G1+G2,3) == 1.0 :
            print(G0+G1+G2)
            GAMMA.append([G0,G1,G2])          
            ADJ.append(adjc)
            
            break
            
t4=time.time()
print(' ***** found exact element and weight function in ', round(t4-t3,3),'s')  
          
Stat_val=[]
## output data
for s in range(len(stat_list)):
    G0,G1,G2=GAMMA[s]
    n0,n1,n2=ADJ[s]
    S=[]
    for t in range(nt):
            mval= G0*var[t][n0]+G1*var[t][n1]+G2*var[t][n2]
            S.append(mval)
    Stat_val.append(S)

t5=time.time()
print(' ***** output ', round(t5-t4,3),'s')  
