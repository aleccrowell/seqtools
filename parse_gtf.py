import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
import pickle

def dictify(l):
    d = {}
    for i in l[:-2].split('; '):
        d[i.split(' ')[0]] = i.split(' ')[1][1:-1]
    return d

def get_keys(dl):
    kl = []
    for i in dl:
        [kl.append(j) for j in list(i.keys())]
    return list(set(kl))

def get_col(dl,k):
    c=[]
    for i in dl:
        try:
            c.append(i[k])
        except:
            c.append('nan')
    return c

def get_chroms(arr):
    d={}
    for i in arr:
        for j in i[1].split(','):
            d[j] = i[0].split('[')[0]
    return d

def check_overlap(row):
    if int(row['Start']) < start < int(row['Stop']):
        return row['ref_gene_id']
    elif int(row['Start']) < stop < int(row['Stop']):
        return row['ref_gene_id']
    elif start < int(row['Start']) < stop:
        return row['ref_gene_id']
    elif start < int(row['Stop']) < stop:
        return row['ref_gene_id']

def get_splice_sites(row):
    nt_starts = merged_exons[merged_exons['transcript_id']==nt]['Start'].tolist()
    nt_stops = merged_exons[merged_exons['transcript_id']==nt]['Stop'].tolist()
    nt_starts.sort()
    nt_stops.sort()



input_fname = 'data/stringtie_merged.gtf'

input_file = open(input_fname, 'r')
skipcount = 0
for i, line in enumerate(input_file):
    if line.startswith('#'):
        skipcount += 1
    else:
        break

in_d = pd.read_csv(input_fname,sep='\t',index_col=False,header=None,skiprows=skipcount)
in_d_2 = pd.read_csv('data/merged.loci',sep='\t',index_col=0,header=None)

int_dict = in_d[8].apply(dictify).values

keys = get_keys(int_dict)

for k in keys:
    in_d[k] = get_col(int_dict,k)

in_d.drop(8,axis=1,inplace=True)
in_d.drop([0,1],axis=1,inplace=True)
in_d.drop(5,axis=1,inplace=True)

in_d.columns = ['Feature_Type','Start','Stop','Strand','Frame']+in_d.columns[5:].tolist()

chroms = get_chroms(in_d_2[[1,3]].values)

chroms = pd.DataFrame.from_dict(chroms,orient='index')
chroms.columns = ['Chrom']

merged = pd.merge(in_d, chroms, left_on='transcript_id', right_index=True)

merged_transcripts = merged[merged['Feature_Type']=='transcript']

######get transcript abundances and filter by actually expressing to limit runtime

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

named.reset_index(inplace=True)
not_novel = named[named['t_name'].str.startswith('ENS')]
novel = named[named['t_name'].str.startswith('MSTR')]
novel.set_index('t_name',inplace=True)
not_novel.set_index('t_name',inplace=True)

novel_real = novel[novel.sum(axis=1)>10]


for t_name in tqdm(novel_real.index.values):
    count = 0
    start = int(merged_transcripts[merged_transcripts['transcript_id']==t_name]['Start'].values[0])
    stop = int(merged_transcripts[merged_transcripts['transcript_id']==t_name]['Stop'].values[0])
    chromosome = merged_transcripts[merged_transcripts['transcript_id']==t_name]['Chrom'].values[0]
    strand = merged_transcripts[merged_transcripts['transcript_id']==t_name]['Strand'].values[0]
    pairs = merged_transcripts[(merged_transcripts['Chrom']==chromosome)&(merged_transcripts['Strand']==strand)&(merged_transcripts['transcript_id'].str.startswith('ENS'))].copy().apply(check_overlap,axis=1).values
    pairs = [i for i in pairs if i != None]
    pairs = list(set(pairs))
    if len(pairs) == 1:
        ind = merged_transcripts[merged_transcripts['transcript_id']==t_name].index.tolist()[0]
        merged_transcripts = merged_transcripts.set_value(ind, 'gene_name', merged_transcripts[merged_transcripts['ref_gene_id']==pairs[0]]['gene_name'].values[0])
        merged_transcripts = merged_transcripts.set_value(ind, 'ref_gene_id', pairs[0])
    elif len(pairs) > 1:
        ind = merged_transcripts[merged_transcripts['transcript_id']==t_name].index.tolist()[0]
        merged_transcripts = merged_transcripts.set_value(ind, 'gene_name', 'multiple')
        merged_transcripts = merged_transcripts.set_value(ind, 'ref_gene_id', 'multiple')


test = merged_transcripts[(merged_transcripts['transcript_id'].isin(novel_real.index.values))&(merged_transcripts['Feature_Type']=='transcript')]




files = [f for f in os.listdir('.') if not os.path.isfile(f)]
files = [f for f in files if f[:2] in ['WT','Wh','UC','Ad']]

temp2 = merged_transcripts.set_index('transcript_id')[['ref_gene_id','gene_name']]
temp2.columns = ['gene_id','gene_name']
temp2 = temp2.replace('nan', np.nan, regex=True)

for f in files:
    temp = pd.read_csv(f+'/t_data.ctab',sep='\t')
    colorder = temp.columns.values
    temp.set_index('t_name',inplace=True)
    temp.update(temp2)
    temp.reset_index(inplace=True)
    temp = temp[colorder]
    temp.to_csv(f+'/t_data.ctab',sep='\t',index=False)


#get overlap dict
overlaps = {}
for t_name in tqdm(novel_real.index.values):
    start = int(merged_transcripts[merged_transcripts['transcript_id']==t_name]['Start'].values[0])
    stop = int(merged_transcripts[merged_transcripts['transcript_id']==t_name]['Stop'].values[0])
    chromosome = merged_transcripts[merged_transcripts['transcript_id']==t_name]['Chrom'].values[0]
    strand = merged_transcripts[merged_transcripts['transcript_id']==t_name]['Strand'].values[0]
    pairs = merged_transcripts[(merged_transcripts['Chrom']==chromosome)&(merged_transcripts['Strand']==strand)&(merged_transcripts['transcript_id'].str.startswith('ENS'))].copy().apply(check_overlap,axis=1).values
    pairs = [i for i in pairs if i != None]
    pairs = list(set(pairs))
    for gene in pairs:
        if gene not in overlaps.keys():
            overlaps[gene] = [t_name]
        else:
            overlaps[gene].append(t_name)

merged_exons = merged[merged['Feature_Type']=='exon']

temp_list = list(set(merged_exons['ref_gene_id'].tolist()))
splice_sites = {}
for i in tqdm(range(len(temp_list))):
    g_name = temp_list.pop()
    current_gene = merged_exons[merged_exons['ref_gene_id']==g_name]
    splice_sites[g_name] = {}
    for t_name in current_gene['transcript_id'].tolist():
        splice_sites[g_name][t_name] = {}
        splice_sites[g_name][t_name]['starts'] = current_gene[current_gene['transcript_id']==t_name]['Start'].tolist()
        splice_sites[g_name][t_name]['stops'] = current_gene[current_gene['transcript_id']==t_name]['Stop'].tolist()

splice_pairs = {}
for g_name in tqdm(overlaps.keys()):
    current_gene = merged_exons[merged_exons['ref_gene_id']==g_name]
    real_ts = current_gene['transcript_id'].tolist()
    for t_name in splice_sites[g_name].keys():
        for nt in overlaps[g_name]:
            nt_starts = merged_exons[merged_exons['transcript_id']==nt]['Start'].tolist()
            nt_stops = merged_exons[merged_exons['transcript_id']==nt]['Stop'].tolist()
            nt_starts.sort()
            nt_stops.sort()
            if len(nt_starts) == len(splice_sites[g_name][t_name]['starts']):
                if 0 == np.sum([i%3 for i in np.subtract(nt_starts[1:],splice_sites[g_name][t_name]['starts'][1:])]) + np.sum([j%3 for j in np.subtract(nt_stops[:-1], splice_sites[g_name][t_name]['stops'][:-1])]):
                    score = np.sum(np.abs([i//3 for i in np.subtract(nt_starts[1:],splice_sites[g_name][t_name]['starts'][1:])])) + np.sum(np.abs([j//3 for j in np.subtract(nt_stops[:-1], splice_sites[g_name][t_name]['stops'][:-1])]))
                    if nt not in splice_pairs.keys():
                        splice_pairs[nt] = [[t_name, score]]
                    elif splice_pairs[nt][0][1] > score:
                        splice_pairs[nt] = [[t_name, score]]
                    elif splice_pairs[nt][0][1] == score:
                        splice_pairs[nt].append([t_name,score])

real_splice_pairs = {}
for key, value in splice_pairs.items():
    if value[0][1] < 3:
        real_splice_pairs[key] = [i[0] for i in value]


pickle.dump( real_splice_pairs, open( "real_splice_pairs.p", "wb" ) )
pickle.dump( splice_pairs, open( "splice_pairs.p", "wb" ) )
pickle.dump( splice_sites, open( "splice_sites.p", "wb" ) )
pickle.dump( overlaps, open( "overlaps.p", "wb" ) )
