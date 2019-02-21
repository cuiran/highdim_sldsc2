import pandas as pd
from subprocess import Popen, PIPE, call

def run():
    phenodf = pd.read_csv('/n/groups/price/ran/high-dim-sldsc/real_pheno_pipeline/ref_files/phenos.csv',delim_whitespace=True)
    n = len(phenodf)
    for i in range(1,3):
        name = phenodf.iloc[i,0]
        ssdir = phenodf.iloc[i,1]
        for j in range(1,23):
                print('running regression leaving out chromosome '+str(j)+' for phenotype '+name)
                cmd = ['sbatch -p short -t 0-12:00 --mem=68000 -o ../output/'+name+'_chr'+str(j)+'_allv2.1.%N_%j.out -e ../output/'+name+'_chr'+str(j)+'_allv2.1.%N_%j.err', 'fitall.sh', ssdir, str(j), name]
                call(' '.join(cmd),shell=True)
    return

if __name__=='__main__':
    run()
