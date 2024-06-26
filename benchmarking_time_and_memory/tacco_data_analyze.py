import scanpy as sc 
import numpy as np 

adata = sc.read_h5ad('tacco_output.h5ad')
cellname=adata.obs_names
print(adata)

cluster1= adata.obsm['Tacco_singleCell_OT']
cluster2= adata.obsm['Tacco_singleCell_withoutmulticenter']

mat1= cluster1.to_numpy()
mat2= cluster2.to_numpy()

print('1',cluster1.shape, type(cluster1))

def find_annot(cluster,header):
    ctname=[]
    d={}
    for i in range(len(cluster)):
        index=np.argmax(cluster[i,:])
        value=header[index]
        ctname.append(index)
        if value in d:
            d[value]+=1
        else:
            d[value]=1
        
    return ctname,d 

ct1,d1= find_annot(mat1,cluster1.columns)
ct2,d2= find_annot(mat2,cluster2.columns)

print(ct1[0:5])

header=cluster1.columns
fw=open('scRNAseq_ctname_ot.dat','w')
for i in range(len(header)):
    fw.write(str(i)+'\t'+header[i]+'\t'+str(d1[header[i]])+'\n')
fw.close()

header=cluster2.columns
fw=open('scRNAseq_ctname_wmc.dat','w')
for i in range(len(header)):
    fw.write(str(i)+'\t'+header[i]+'\t'+str(d2[header[i]])+'\n')
fw.close()

fw=open('tacoo_derived_cluster.dat','w')
fw.write('barcode,tacco,cluster\n')
for i in range(len(cellname)):
    fw.write(cellname[i]+','+str(ct1[i])+','+str(ct2[i])+'\n')
fw.close()
