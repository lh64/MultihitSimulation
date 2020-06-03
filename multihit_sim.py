import numpy as np
import matplotlib.pyplot as plt
import re
import os

filename = os.path.join(os.path.curdir, 'SIM_INPUT', 'gene_length_mut_50strains.csv')
suffix = re.match(".+mut(.*)\.csv",filename).group(1)
data = np.genfromtxt(filename, dtype=None, delimiter=',', names=True, encoding='UTF-8')

print(data.dtype.names)

HO_nonsyn = [(d['Length'], d['Nonsyn']) for d in data 
             if d['Orientation'] == 'HO' and d['Length'] > 200]
CD_nonsyn = [(d['Length'], d['Nonsyn']) for d in data 
             if d['Orientation'] == 'CD' and d['Length'] > 200]

## HO and CD genes have same length distributions 
'''
bins = np.arange(0,2001,100)
plt.hist([ho[0] for ho in HO_nonsyn], bins=bins, label='HO', 
         weights=[1/len(HO_nonsyn)]*len(HO_nonsyn))
plt.hist([cd[0] for cd in CD_nonsyn], bins=bins, label='CD', fc='None', ec='k', lw=2,
         weights=[1/len(CD_nonsyn)]*len(CD_nonsyn))
plt.show()
'''

np.random.seed(0)

# percentage of variable sites
pct_var = np.arange(0.1, 1.1, 0.1)

for pct in pct_var:
    
    niter = 10000
    
    sum_hsmut_sites_HO = [0]*niter
    sum_hsmut_sites_CD = [0]*niter
    sum_hsmut_genes_HO = [0]*niter
    sum_hsmut_genes_CD = [0]*niter
    sum_parmut_HO  = [0]*niter
    sum_parmut_CD  = [0]*niter
    
    for n in range(niter):
        print(n)
        sum_hsmut_sites = [0,0] # [HO,CD]
        sum_hsmut_genes = [0,0]
        sum_parmut = [0,0]
        for i,type in enumerate([HO_nonsyn, CD_nonsyn]):
            for gene in type:
                nsites = int(round(gene[0]*pct))
                x = np.random.randint(nsites, size=gene[1])
                pos, hits = np.unique(x, return_counts=True)
                multihits = [h for h in hits if h > 1]
                n_multihits = len(multihits)
    
                if n_multihits > 0:
                    sum_hsmut_sites[i] += n_multihits
                    sum_hsmut_genes[i] += 1
                    # Check for parallel mutations
                    for mh in multihits:
                        # 20 amino acids to choose from (there are 21 total)
                        xx = np.random.randint(20, size=mh) 
                        aa, nmut = np.unique(xx, return_counts=True)
                        parallel = [nm for nm in nmut if nm > 1]
                        n_parallel = len(parallel)
                        if n_parallel > 0:
                            sum_parmut[i] += n_parallel
                
        sum_hsmut_sites_HO[n] = sum_hsmut_sites[0]
        sum_hsmut_sites_CD[n] = sum_hsmut_sites[1]
        sum_hsmut_genes_HO[n] = sum_hsmut_genes[0]
        sum_hsmut_genes_CD[n] = sum_hsmut_genes[1]
        sum_parmut_HO[n]  = sum_parmut[0]
        sum_parmut_CD[n]  = sum_parmut[1]
        
    outfile = os.path.join(os.path.curdir, 'SIM_OUTPUT', 
                           'hsmut_sites%s_%d.csv' % (suffix, pct*100))
    with open(outfile, 'w') as f:
        f.write("HO,CD\n")
        for n in range(niter):
            f.write("%d,%d\n" % (sum_hsmut_sites_HO[n],sum_hsmut_sites_CD[n]))
    
    outfile = os.path.join(os.path.curdir, 'SIM_OUTPUT', 
                           'hsmut_genes%s_%d.csv' % (suffix, pct*100))
    with open(outfile, 'w') as f:
        f.write("HO,CD\n")
        for n in range(niter):
            f.write("%d,%d\n" % (sum_hsmut_genes_HO[n],sum_hsmut_genes_CD[n]))
    
    outfile = os.path.join(os.path.curdir, 'SIM_OUTPUT', 
                           'parmut%s_%d.csv' % (suffix, pct*100))
    with open(outfile, 'w') as f:
        f.write("HO,CD\n")
        for n in range(niter):
            f.write("%d,%d\n" % (sum_parmut_HO[n],sum_parmut_CD[n]))
