# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 12:08:41 2021

@author: nicol
"""

import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', size=14)
rc('font', family='Arial')
# rc('font', serif='Times')
rc('pdf', fonttype=42)
# rc('text', usetex=True)

folder = os.path.join('Y:',os.sep,'Nicola_Gritti','analysis_code','scRNAseq_Gastruloids','new_codes','results','integration','pijuan_anlas')

umap = pd.read_csv(os.path.join(folder,'integration_umap.csv'))
umap['cellname'] = umap.index
umap = umap.reset_index()
meta_anlas = pd.read_csv(os.path.join(folder,'meta_anlas.csv'))
meta_pijuan = pd.read_csv(os.path.join(folder,'meta_pijuan.csv'))
meta = pd.concat([meta_pijuan,meta_anlas], ignore_index=True)
meta = meta[['celltype.general','batch.ident','celltype.anlas','celltype.pijuan','stage','merge.ident']]

meta = meta.rename(columns={'celltype.general':'celltype',
                            'celltype.anlas':'celltype_anlas',
                            'celltype.pijuan':'celltype_pijuan',
                            'merge.ident':'stage_anlas',
                            'stage':'stage_pijuan',
                            'batch.ident':'batch_ident'})
reps = []
for i in meta.batch_ident:
    if 'rep1' in i:
        reps.append('rep1')
    elif 'rep2' in i:
        reps.append('rep2')
    else:
        reps.append(np.nan)
meta['replicate'] = reps

meta = meta.replace(np.nan, 'None', regex=True)
stage = []
for i in range(len(meta.celltype)):
    sa = meta.stage_anlas[i]
    sp = meta.stage_pijuan[i]
    if sa=='None':
        stage.append(sp)
    else:
        stage.append(sa)
meta['stage'] = stage

df = pd.concat([umap,meta], axis=1)

fig, ax = plt.subplots(figsize=(6,6))
fig.subplots_adjust(left=0.15, right=0.99, top=0.99, bottom=0.15)
sns.scatterplot(data=df, x='V1', y='V2', hue='celltype', 
                s=2,
                hue_order=['Epiblast','PGC','ExE','Endoderm','Ectoderm','Mesoderm','Primitive Streak',
                           'G Pluripotent','G Primed pluripotent','G Early diff','G Endodermal','G Neural','G Mesoderm','G Mesendodermal'],
                linewidth=0, alpha=.4, ax=ax)
plt.legend(loc='upper right', frameon=False)
