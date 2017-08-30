import pandas as pd
import numpy as np
import os

def qnorm(df):
    ref = pd.concat([df[col].sort_values().reset_index(drop=True) for col in df], axis=1, ignore_index=True).mean(axis=1).values
    for i in range(0,len(df.columns)):
        df = df.sort_values(df.columns[i])
        df[df.columns[i]] = ref
    return df.sort_index()

files = [f for f in os.listdir('.') if not os.path.isfile(f)]
files = [f for f in files if f[:2]=='WT']

temp = pd.read_csv(files[0]+'/t_data.ctab',sep='\t',index_col=0)['FPKM'].to_frame()
temp.columns = ['ZT'+files[0].split('_')[0].split('WT')[1]]
for f in files[1:]:
    temp2 = pd.read_csv(f+'/t_data.ctab',sep='\t',index_col=0)['FPKM'].to_frame()
    temp2.columns = ['ZT'+f.split('_')[0].split('WT')[1]]
    temp = temp.merge(temp2, right_index=True, left_index=True, how='outer')
    print('ZT'+f.split('_')[0].split('WT')[1]+' OK')
    print(len(temp))

names = pd.read_csv(files[0]+'/t_data.ctab',sep='\t',index_col=None)[['t_id','t_name']]
names.set_index('t_id',inplace=True)

temp = temp.merge(names, right_index=True, left_index=True, how='outer')
named = temp.set_index('t_name')

qnormed = qnorm(named)
qnormed.index.names = ['#']
qnormed = qnormed[qnormed.sum(axis=1)>10]
qnormed = np.log(qnormed+1)
qnormed.to_csv('WT_transcript_for_ejtk.txt',sep='\t')
