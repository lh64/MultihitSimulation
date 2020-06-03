import numpy as np
import matplotlib.pyplot as plt
import os

SUFFIX = ['_5strains', '_50strains']
EC = ['black', 'black']
FC = ['0.75','None']
LW = [1,2]
LABEL = ['5 strains', '50 strains']
PCT_VAR = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

def make_plots(SUFFIX,EC,FC,LW,LABEL):
    
    for pct in PCT_VAR:
            
        print(pct)
    
        for suffix,ec,fc,lw,label in zip(SUFFIX,EC,FC,LW,LABEL):
            
            print(suffix)
 
            ##### MULTI-HIT MUTATIONS #####
    
            filename = os.path.join(os.path.curdir, 'SIM_INPUT', 
                                    'gene_length_mut%s.csv' % suffix)
            data_exp = np.genfromtxt(filename, dtype=None, delimiter=',', names=True, 
                                     encoding='UTF-8')
            
            HO_mut_exp = [(d['HSmut'], d['Convergent'], d['Parallel'], d['Coincidental']) 
                      for d in data_exp if d['Orientation'] == 'HO' and d['Length'] > 200]
            CD_mut_exp = [(d['HSmut'], d['Convergent'], d['Parallel'], d['Coincidental']) 
                      for d in data_exp if d['Orientation'] == 'CD' and d['Length'] > 200]
            
            n_ho_genes = len(HO_mut_exp)
            n_cd_genes = len(CD_mut_exp)
                    
            ho_hsmut_exp = [x[0] for x in HO_mut_exp if x[0]>0]
            cd_hsmut_exp = [x[0] for x in CD_mut_exp if x[0]>0]
            
            ratio_sites_exp = sum(ho_hsmut_exp)/sum(cd_hsmut_exp)
            ratio_genes_exp = len(ho_hsmut_exp)/len(cd_hsmut_exp)
            
            filename = os.path.join(os.path.curdir, 'SIM_OUTPUT', 
                                    'hsmut_sites%s_%d.csv' % (suffix,pct))
            data_sim_sites = np.genfromtxt(filename, dtype=None, delimiter=',', names=True, 
                                           encoding='UTF-8')
            filename = os.path.join(os.path.curdir, 'SIM_OUTPUT', 
                                    'hsmut_genes%s_%d.csv' % (suffix,pct))
            data_sim_genes = np.genfromtxt(filename, dtype=None, delimiter=',', names=True, 
                                           encoding='UTF-8')
            
            ho_hsmut_sites_sim = [d[0] for d in data_sim_sites]
            cd_hsmut_sites_sim = [d[1] for d in data_sim_sites]
            ratio_sites_sim = [d[0]/d[1] for d in data_sim_sites]
            
            ho_hsmut_genes_sim = [d[0] for d in data_sim_genes]
            cd_hsmut_genes_sim = [d[1] for d in data_sim_genes]
            ratio_genes_sim = [d[0]/d[1] for d in data_sim_genes]
            
            def get_pval(sim_data, exp_val):
                p_geq = len([r for r in sim_data if r >= exp_val]) / \
                    len(sim_data)
                p_leq = len([r for r in sim_data if r <= exp_val]) / \
                    len(sim_data)
                return 2*min(p_geq, p_leq) #two-tailed p-value
            
            plt.figure('R_N_%d' % pct)
            dx = 0.025 #0.1
            bins = np.arange(0, 0.8+dx, dx)
            plt.hist(ratio_sites_sim, bins=bins, density=False, 
                     weights=[1/len(ratio_sites_sim)]*len(ratio_sites_sim), 
                     color='white', ec=ec, fc=fc, lw=lw, label=label)
            p_value = get_pval(ratio_sites_sim, ratio_sites_exp)
            plt.plot([ratio_sites_exp], [0.3], '*', ms=12, mec=ec, mfc=fc, 
                     label='Observed (p = %.2f)' % p_value)
            plt.xlabel(r'$R_N$')
            plt.ylabel('frequency')
            plt.legend(loc=0)
    
            plt.figure('R_G_%d' % pct)
            dx = 0.025 #0.1
            bins = np.arange(0, 0.8+dx, dx)
            plt.hist(ratio_genes_sim, bins=bins, density=False, 
                     weights=[1/len(ratio_genes_sim)]*len(ratio_genes_sim), 
                     color='white', ec=ec, fc=fc, lw=lw, label=label)
            p_value = get_pval(ratio_genes_sim, ratio_genes_exp)
            plt.plot([ratio_genes_exp], [0.3], '*', ms=12, mec=ec, mfc=fc, 
                     label='Observed (p = %.2f)' % p_value)
            plt.xlabel(r'$R_G$')
            plt.ylabel('frequency')
            plt.legend(loc=0)
            
            def get_bins(x, y):
                dmin = min(min(x),min(y))
                dmax = max(max(x),max(y))
                order = 0
                n = 2
                while dmax > 2*10**(order+(order+1)*n):
                    order += 1
                step = 10**order
                if dmin % 10*step != 0:
                    dmin = int(round((dmin-5*step)/(10*step))*(10*step))
                if dmax % 10*step != 0:
                    dmax = int(round((dmax+5*step)/(10*step))*(10*step))
                bins = range(dmin,dmax,step)
                return bins
            
            plt.figure('hsmut_sites%s_%d' % (suffix,pct))
            bins = get_bins(ho_hsmut_sites_sim, cd_hsmut_sites_sim)
            plt.hist(ho_hsmut_sites_sim, bins=bins, color='white', fc='0.75', 
                     ec='red', lw=1, label='simulated HO')
            plt.hist(cd_hsmut_sites_sim, bins=bins, color='white', fc='None', 
                     ec='blue', lw=2, label='simulated CD')
            p_value_ho = get_pval(ho_hsmut_sites_sim, sum(ho_hsmut_exp))
            p_value_cd = get_pval(cd_hsmut_sites_sim, sum(cd_hsmut_exp))
            plt.annotate('Observed HO = {:,} (p = {:.4g})'.format(sum(ho_hsmut_exp), p_value_ho), 
                         xy=(0.4,0.75), xycoords='axes fraction', weight='bold', color='red')
            plt.annotate('Observed CD = {:,} (p = {:.4g})'.format(sum(cd_hsmut_exp), p_value_cd), 
                         xy=(0.4,0.69), xycoords='axes fraction', weight='bold', color='blue')
            plt.xlabel('# multi-hit sites') 
            plt.ylabel('number')
            plt.legend(loc=0)      
            
            outfile = os.path.join(os.path.curdir, 'PLOTS', 
                                   'hsmut_sites%s_%d.pdf' % (suffix,pct))
            plt.savefig(outfile, format='pdf') 
    
            plt.figure('hsmut_genes%s_%d' % (suffix,pct))
            bins = get_bins(ho_hsmut_genes_sim, cd_hsmut_genes_sim)
            plt.hist(ho_hsmut_genes_sim, bins=bins, color='white', fc='0.75', 
                     ec='red', lw=1, label='simulated HO')
            plt.hist(cd_hsmut_genes_sim, bins=bins, color='white', fc='None', 
                     ec='blue', lw=2, label='simulated CD')
            p_value_ho = get_pval(ho_hsmut_genes_sim, len(ho_hsmut_exp))
            p_value_cd = get_pval(cd_hsmut_genes_sim, len(cd_hsmut_exp))
            plt.annotate('Observed HO = {:,} (p = {:.4g})'.format(len(ho_hsmut_exp), p_value_ho), 
                         xy=(0.4,0.75), xycoords='axes fraction', weight='bold', color='red')
            plt.annotate('Observed CD = {:,} (p = {:.4g})'.format(len(cd_hsmut_exp), p_value_cd), 
                         xy=(0.4,0.69), 
                         xycoords='axes fraction', weight='bold', color='blue')
            plt.xlabel('# multi-hit genes') 
            plt.ylabel('number')
            plt.legend(loc=0)
            
            outfile = os.path.join(os.path.curdir, 'PLOTS', 
                                   'hsmut_genes%s_%d.pdf' % (suffix,pct))
            plt.savefig(outfile, format='pdf')   
            
            ##### PARALLEL MUTATIONS #####
            
            ho_par_exp = [x[2] for x in HO_mut_exp if x[2]>0]
            cd_par_exp = [x[2] for x in CD_mut_exp if x[2]>0]
            
            ratio_sites_par_exp = sum(ho_par_exp)/sum(cd_par_exp)
            ratio_genes_par_exp = len(ho_par_exp)/len(cd_par_exp)
            
            filename = os.path.join(os.path.curdir, 'SIM_OUTPUT', 
                                    'parmut%s_%d.csv' % (suffix,pct))
            data_sim_par = np.genfromtxt(filename, dtype=None, delimiter=',', names=True, 
                                         encoding='UTF-8')
            
            sum_par_HO = [d[0] for d in data_sim_par]
            sum_par_CD = [d[1] for d in data_sim_par]
            
            plt.figure('parmut%s_%d' % (suffix,pct))
            bins = get_bins(sum_par_HO, sum_par_CD)
            plt.hist(sum_par_HO, bins=bins, color='white', fc='0.75', 
                     ec='red', lw=1, label='simulated HO')
            plt.hist(sum_par_CD, bins=bins, color='white', fc='None', 
                     ec='blue', lw=2, label='simulated CD')
            p_value_ho = get_pval(sum_par_HO, sum(ho_par_exp))
            p_value_cd = get_pval(sum_par_CD, sum(cd_par_exp))
            plt.annotate('Observed HO = {:,} (p = {:.4g})'.format(sum(ho_par_exp), p_value_ho), 
                         xy=(0.4,0.75), xycoords='axes fraction', weight='bold', color='red')
            plt.annotate('Observed CD = {:,} (p = {:.4g})'.format(sum(cd_par_exp), p_value_cd), 
                         xy=(0.4,0.69), xycoords='axes fraction', weight='bold', color='blue')
            plt.xlabel('# parallel mutations') 
            plt.ylabel('number')
            plt.legend(loc=0)
            
            outfile = os.path.join(os.path.curdir, 'PLOTS', 
                                   'parmut%s_%d.pdf' % (suffix,pct))
            plt.savefig(outfile, format='pdf')
        
        plt.figure('R_N_%d' % pct)
        ymax = round(plt.ylim()[1]*10+0.5)/10
        plt.ylim(ymax=ymax) #0.5
        outfile = os.path.join(os.path.curdir, 'PLOTS', 'multihit_RN_%d.pdf' % pct)
        plt.savefig(outfile, format='pdf')
        
        plt.figure('R_G_%d' % pct)
        ymax = round(plt.ylim()[1]*10+0.5)/10
        plt.ylim(ymax=ymax) #0.5
        outfile = os.path.join(os.path.curdir, 'PLOTS', 'multihit_RG_%d.pdf' % pct)
        plt.savefig(outfile, format='pdf')
        
#         plt.show()
    
make_plots(SUFFIX, EC, FC, LW, LABEL)
