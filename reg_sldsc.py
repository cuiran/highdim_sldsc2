import pandas as pd
import numpy as np
from sklearn import linear_model
import argparse
import pdb
from sklearn.metrics import r2_score
import h5py
from keras.models import Sequential
from keras.models import Model
from keras.layers.core import Dense
from keras import optimizers
from keras import regularizers
import shrink
import sys

def get_data_singlefile(args):
    ld_file = args.ld
    sumstat_file = args.sumstat
    weightld_file = args.weight_ld
    out_chr = args.leave_chrom
    outfile_prefix = args.outfile_prefix

    print('reading in data')
    ld_df = pd.read_csv(ld_file,delim_whitespace=True)
    ss_df = pd.read_csv(sumstat_file,delim_whitespace=True)
    wld_df = pd.read_csv(weightld_file,delim_whitespace=True)
    if int(out_chr) == 0:
        if args.reg_weights:
            # this means regression weights is given
            # reg_weights is approx reciprocal of conditional variance function
            print('read in regression weights')
            w = np.array(pd.read_csv(args.reg_weights,delim_whitespace=True)['WEIGHT'])
            w = 1.0/w

    N = np.mean(ss_df['N'])
    ld_wld_merged = pd.merge(ld_df,wld_df[['SNP','L2']],on=['SNP']) 
    print('merging data')
    all_merged = pd.merge(ld_wld_merged,ss_df[['SNP','CHISQ']],on=['SNP'])
    chisq_max = 80.0 #max(0.001*max(ss_df['N']),80)
    all_merged = all_merged[all_merged['CHISQ']<chisq_max]
    print('after merging, num SNPs remain:', len(all_merged))
    y = all_merged['CHISQ']
    wld = all_merged['L2']
    merged_ld = all_merged.drop(['CHR','SNP','BP','CHISQ','L2'],axis=1)
    X = np.array(merged_ld.iloc[:,:])
    if not args.reg_weights:
        print('computing regression weights')
        #M = get_M(args.annot_prefix,'.l2.M')
        SNP_info = all_merged[['CHR','SNP','BP']].copy()
        w = compute_weights(X,y,wld,outfile_prefix,SNP_info)
        if args.just_weights:
            sys.exit()
    all_merged['L2'] = w

    out_chr = int(out_chr)
    if out_chr == 0:
        train_df = all_merged.copy()
    else:
        train_df = all_merged[all_merged['CHR']!=out_chr]
        test_df = all_merged[all_merged['CHR']==out_chr]

    print('spliting training and testing data')
    y_train = np.array(train_df['CHISQ'])
    w_train = np.array(train_df['L2'])
    train_ld = train_df.drop(['CHR','SNP','BP','CHISQ','L2'],axis=1)
    X_train = np.array(train_ld.iloc[:,:])

    if out_chr == 0:
        y_test = 0
        X_test = 0
        w_test = 0
    else:
        y_test = np.array(test_df['CHISQ'])
        w_test = np.array(test_df['L2'])
        test_ld = test_df.drop(['CHR','SNP','BP','CHISQ','L2'],axis=1)
        X_test = np.array(test_ld.iloc[:,:])
    return X_train,X_test,y_train,y_test,w_train,w_test,N


def get_M(mprefix,msuffix):
    # assuming there's only one M file per chromosome
    # output an array of length equal to num annots
    # with each element the sum of M files across all chromosomes
    fnames = [mprefix+str(i)+msuffix for i in range(1,23)]
    mdfs = [pd.read_csv(f,delim_whitespace=True,header=None) for f in fnames]
    mdf = pd.concat(mdfs)
    m = np.array(mdf.sum())
    return m

def compute_weights(ld,ss,wld,outfile_prefix,SNP_info_df):
    # return regression weights, approximately the conditional variance
    wld = np.array(np.fmax(wld,1.0))
    #M = ld.shape[0]
    #pdb.set_trace()
    #sum_ld = np.sum(ld)
    #sum_ss = np.sum(ss)
    #l = sum_ld/M
    #s = sum_ss/M
    #Ntau_hat = (s-1)/l
    #sum_ld_col = np.sum(ld,axis=1)
    #sum_ld_col = np.fmax(sum_ld_col,1.0)
    #het_w = 2*((Ntau_hat*sum_ld_col +1)**2)
    #het_w = np.fmax(het_w,1.0)
    #w = np.multiply(het_w,wld)
    w = wld
    #SNP_info_df['WEIGHT'] = w
    #SNP_info_df.to_csv(outfile_prefix+'_reg_weights.csv.tmp',sep='\t',index=False)
    #print('weights stored at '+outfile_prefix+'_reg_weights.csv.tmp')
    return w

def compute_weights2(ld,ss,wld,outfile_prefix,M,N):
    # implemented according to Omer's code, but the result looks more different from s-ldsc weights than mine does
    sum_ld_col = np.sum(ld,axis=1)
    hsq = ((ss.mean() - 1) * M.sum()/N) / np.mean(sum_ld_col)
    hsq = np.clip(hsq, 0, 1)
    sum_ld_col = np.fmax(sum_ld_col,1.0)
    wld = np.array(np.fmax(wld,1.0))
    c = hsq * N / M.sum()
    intercept_temp = 1.0
    het_w = 1.0 / (2 * (intercept_temp + c*sum_ld_col)**2)
    oc_w = 1.0 / wld
    ldsc_w = np.sqrt(het_w * oc_w)
    ldsc_w /= ldsc_w.sum()
    return ldsc_w
    

def fit_lasso(X,y):
    print('performing Lasso regression')
    model = linear_model.LassoCV(max_iter = 1000000, fit_intercept=False)
    model.fit(X,y)
    return model

def fit_ols(X,y):
    print('performing OLS regression')
    model = linear_model.LinearRegression(fit_intercept=False)
    model.fit(X,y)
    return model

def fit_elnet(X,y):
    print('performing Elastic Net regression')
    model = linear_model.ElasticNetCV(max_iter = 1000000, l1_ratio = [.001, .1, .5, .7, .9, .95, .99, 1],fit_intercept=False)
    model.fit(X,y)
    return model

def fit_ridge(X,y,lambdas):
    print('performing Ridge regression')
    model = linear_model.RidgeCV(alphas = lambdas,fit_intercept=False)
    model.fit(X,y)
    return model

def keras_lasso(X,y,rp,lr_hat):
    # rp is regularization parameter, here we use the alpha generated by sklearn LassoCV
    print('performing Lasso regression using keras')
    model = Sequential()
    p = X.shape[1]
    model.add(Dense(1,input_dim=p,use_bias=False))
    sgd = optimizers.SGD(lr=lr_hat)
    model.compile(loss='mse',optimizer=sgd)
    batch_size = 1000
    steps = X.shape[0]//batch_size
    l1_shrink = shrink.L1_update(model.trainable_weights[:1],lr=lr_hat,regularizer=rp)
    new_X = (X-np.mean(X,axis=0))/np.std(X,axis=0)
    new_y = (y-np.mean(y))/np.std(y)
    model.fit(new_X,new_y,batch_size=batch_size,epochs=100,
                        verbose=1,callbacks=[l1_shrink])
    coefs = model.get_weights()[0]
    return coefs.reshape((coefs.shape[0],))

def est_lr(N,num_SNPs,num_features):
    tau_hat = N/(num_SNPs*num_features)
    lr_hat = 10**np.round(np.log10(tau_hat/5))
    return lr_hat

def recover(learned_coef,N,meanX,meany):
    true_coef = learned_coef/N 
    c = meany - meanX.dot(learned_coef) 
    return true_coef,c 

def recover2(learned_coef,N,meanX,meany,std_wcX,std_wcy):
    true_coef = np.multiply(std_wcy/std_wcX,learned_coef)
    c = meany - meanX.dot(true_coef)
    tau = true_coef/N
    return tau,c

def generator(X,y,m):
    # X is regressor, y is target, m is batch size
    M = X.shape[0]
    n = M//m #number of complete batches in 1 iteration of X
    while True:
        for i in range(1,n+2): # the last batch would be the combination of the tail of the data and the head
            if i<n+1:
                bx = X[m*(i-1):m*i,:] #batch X
                by = y[m*(i-1):m*i] #batch y
            else:
                n_top = m - M%m #number of lines needed from top
                bx = np.concatenate((X[m*(i-1):,:],X[0:n_top]),axis=0)
                by = np.concatenate((y[m*(i-1):],y[0:n_top]),axis=0)
            yield bx,by

def stdize1(X,y,w):
    # standardize and weight training data
    # return regX, regy, meanX, meany, sqrt_inv_w
    meanX = np.mean(X,axis=0)
    meany = np.mean(y)
    cX = X - meanX
    cy = y - meany
    sqrt_inv_w = np.sqrt(1/w)
    wcX = sqrt_inv_w[:,None]*cX
    wcy = sqrt_inv_w[:,None]*cy[:,None]
    regX = wcX
    regy = wcy.ravel()
    ########### center wcX and wcy again with sum of sqrt_inv_w ########
    regX = regX/np.sum(sqrt_inv_w)
    regy = regy/np.sum(sqrt_inv_w)
    return regX,regy,meanX,meany,sqrt_inv_w

def stdize2(X,y,w):
    # center-weight-scale
    # return regX, regy, meanX, meany, sqrt_inv_w, std_wcX, std_wcy
    meanX = np.mean(X,axis=0)
    meany = np.mean(y)
    cX = X - meanX
    cy = y - meany
    sqrt_inv_w = np.sqrt(1/w)
    wcX = sqrt_inv_w[:,None]*cX
    wcy = sqrt_inv_w[:,None]*cy[:,None]
    std_wcX = np.std(wcX,axis=0)
    std_wcy = np.std(wcy)
    regX = wcX/std_wcX
    regy = wcy.ravel()/std_wcy
    return regX, regy, meanX, meany, sqrt_inv_w, std_wcX, std_wcy


def run_reg(args):
    outfile_prefix = args.outfile_prefix
    X_train,X_test,y_train,y_test,w_train,w_test,N = get_data_singlefile(args)
    #regX,regy,meanX,meany,sqrt_inv_w = stdize1(X_train,y_train,w_train)
    regX,regy,meanX,meany,sqrt_inv_w,std_wcX,std_wcy = stdize2(X_train,y_train,w_train)
    for method in args.methods:
        if method=='OLS':
            model = fit_ols(regX,regy)
            coef = model.coef_
        elif method=='Lasso':
            model = fit_lasso(regX,regy)
            alphadf = pd.DataFrame(data=[model.alpha_],columns=['ALPHA'])
            alphadf.to_csv(args.outfile_prefix+'_Lasso_alpha',index=False,sep='\t')
            coef = model.coef_
        elif method=='Elnet':
            model = fit_elnet(regX,regy)
            alphadf = pd.DataFrame(data=[model.alpha_],columns=['ALPHA'])
            alphadf.to_csv(args.outfile_prefix+'Elnet_alpha',index=False,sep='\t')
            l1ratiodf = pd.DataFrame(data=[model.l1_ratio_],columns=['L1_Ratio'])
            l1ratiodf.to_csv(args.outfile_prefix+'Elnet_l1_ratio',index=False,sep='\t')
            coef = model.coef_
        elif method=='Ridge':
            mean_XX = get_meanXX(regX)
            print('computing candidate regularization parameters for Ridge regression')
            ridge_lambdas = np.logspace(np.log10(mean_XX*1e-8),np.log10(mean_XX*1e2),num=100)
            model = fit_ridge(regX,regy,ridge_lambdas)
            alphadf = pd.DataFrame(data=[model.alpha_],columns=['ALPHA'])
            alphadf.to_csv(args.outfile_prefix+'_Ridge_alpha',index=False,sep='\t')
            coef = model.coef_
        elif method=='kerasLasso':
            rp = get_alpha(args.alpha_file)
            lr_hat = est_lr(N,X_train.shape[0],X_train.shape[1])
            coef = keras_lasso(regX,regy,rp,lr_hat)
        elif method=='mLasso':
            # mLasso means we fit each annotation marginally first, threshold by significance, then perform joint Lasso fit
            print('marginal+Lasso has not been implemented yet')
        #tau,c = results_to_files1(coef,N,meanX,meany,outfile_prefix,method)
        tau,c = results_to_files2(coef,N,meanX,meany,std_wcX,std_wcy,outfile_prefix,method)
        if args.leave_chrom != 0:
            compute_losses(tau,c,X_train,y_train,X_test,y_test,w_test,N,args.outfile_prefix+'_'+method)
            compute_varbeta(tau,args.leave_chrom,args.annot_prefix,args.outfile_prefix+'_'+method,args.rsid,args.m550_rsid)
    return

def get_alpha(f):
    df = pd.read_csv(f,delim_whitespace=True)
    return df['ALPHA'][0]

def get_meanXX(regX):
    diag = list()
    for i in range(regX.shape[1]):
        x = regX[:,i]
        diag.append(x.dot(x))
    return np.mean(diag)

def results_to_files2(learned_coef,N,meanX,meany,std_wcX,std_wcy,outfile_prefix,method):
    print('saving results')
    tau,c = recover2(learned_coef,N,meanX,meany,std_wcX,std_wcy)
    taudf = pd.DataFrame(data=tau,columns=['TAU'])
    taudf.to_csv(outfile_prefix+'_'+method+'_coef',index=False,sep='\t')
    cdf = pd.DataFrame(data=[c],columns=['INTERCEPT'])
    cdf.to_csv(outfile_prefix+'_'+method+'_intercept',index=False,sep='\t')
    return tau,c

def results_to_files1(learned_coef,N,meanX,meany,outfile_prefix,method):
    print('saving results')
    tau,c = recover(learned_coef,N,meanX,meany)
    taudf = pd.DataFrame(data=tau,columns=['TAU'])
    taudf.to_csv(outfile_prefix+'_'+method+'_coef',index=False,sep='\t')
    cdf = pd.DataFrame(data=[c],columns=['INTERCEPT'])
    cdf.to_csv(outfile_prefix+'_'+method+'_intercept',index=False,sep='\t')
    return tau,c

def compute_varbeta(tau,chrom,annot_prefix,outfile_prefix,rsid_prefix,m550rsid_prefix):
    # annotation prefix + .annot.h2
    # rsid prefix + i
    chrom = str(chrom)
    print('computing per-SNP heritability on chromosome '+chrom)
    annot_file = annot_prefix+chrom+'.annot.gz'
    rsfile = rsid_prefix+chrom
    anndf = pd.read_csv(annot_file,delim_whitespace=True,header=None)
    chr_annot = np.array(anndf.iloc[:,:])
    chr_varbeta = chr_annot.dot(tau)
    rsdf = pd.read_csv(rsfile,delim_whitespace=True)
    rsdf['VARBETA'] = chr_varbeta
    rsdf.to_csv(outfile_prefix+'_all_varbeta',sep='\t',index=False)
    all_h2 = np.sum(chr_varbeta)
    allh2df = pd.DataFrame(data=[all_h2],columns=['H2'])
    allh2df.to_csv(outfile_prefix+'_all_h2',sep='\t',index=False)
    m550rs = pd.read_csv(m550rsid_prefix+chrom,delim_whitespace=True)
    m550rs_snp = m550rs[['SNP']]
    merged = pd.merge(rsdf,m550rs_snp,on=['SNP'])
    merged.to_csv(outfile_prefix+'_common_varbeta',sep='\t',index=False)
    h2 = np.sum(list(merged['VARBETA']))
    h2df = pd.DataFrame(data=[h2],columns=['H2'])
    h2df.to_csv(outfile_prefix+'_common_h2',sep='\t',index=False)
    return 


def compute_losses(tau,c,X_train,y_train,X_test,y_test,w_test,N,outfile_prefix):
    ypred = N*X_test.dot(tau)+c
    preddf = pd.DataFrame(ypred,columns=['PRED_CHISQ'])
    preddf.to_csv(outfile_prefix+'_predchisq',index=False,sep='\t')
    mean_est = np.mean(y_train)
        
    inv_w = 1/w_test
    persnp_wsse = inv_w.dot((ypred-y_test)**2)/np.sum(inv_w)
    persnp = pd.DataFrame(data=[persnp_wsse],columns=['pwsse'])
    persnp.to_csv(outfile_prefix+'_pwsse',index=False,sep='\t')
    worst_persnp = inv_w.dot((mean_est-y_test)**2)/np.sum(inv_w)
    worst_psdf = pd.DataFrame(data=[worst_persnp],columns=['pwsse'])
    worst_psdf.to_csv(outfile_prefix+'_mean_est_pwsse',index=False,sep='\t')
    return
   
def compute_best_pwsse(args):
    true_tau_file = args.true_tau_file
    test_file = args.test_file
    outfile_prefix = args.outfile_prefix
    N = 47360

    true_tau = np.array(pd.read_csv(true_tau_file,delim_whitespace=True)['TAU'])
    test_df = pd.read_csv(test_file,delim_whitespace=True)
    X_test = test_df.drop(['SNP','CHR','L2','CHISQ'],axis=1)
    X_test = np.array(X_test.iloc[:,:])
    w_test = np.array(test_df['L2'])
    y_test = np.array(test_df['CHISQ'])
    ypred = N*X_test.dot(true_tau)+1

    inv_w = 1/w_test 
    persnp_wsse = inv_w.dot((ypred-y_test)**2)/np.sum(inv_w)
    persnp = pd.DataFrame(data=[persnp_wsse],columns=['pwsse'])
    persnp.to_csv(outfile_prefix+'_best_pwsse',index=False,sep='\t')
    return
    
def compute_varbeta_loss(args):
    annot_file = args.annot_prefix
    coef_file = args.coef_file
    true_varbeta_file = args.true_varbeta_file
    outfile_prefix = args.outfile_prefix

    annot = h5py.File(annot_file,'r')['dataset'][:]
    coef = np.array(pd.read_csv(coef_file,delim_whitespace=True)['TAU'])
    pred_varbeta = annot.dot(coef)
    varbeta = np.array(pd.read_csv(true_varbeta_file,delim_whitespace=True)['VARBETA'])
    sse = np.sum((pred_varbeta - varbeta)**2)

    ssedf = pd.DataFrame(data=[sse],columns=['varbeta_sse'])
    ssedf.to_csv(outfile_prefix+'_varbeta_sse',index=False,sep='\t')
    return

def correct_pwsse(args):
    coef_file = args.coef_file
    test_file = args.test_file
    N = 47360
    outfile_prefix = args.outfile_prefix

    tau = np.array(pd.read_csv(coef_file,delim_whitespace=True)['TAU'])
    test_df = pd.read_csv(test_file,delim_whitespace=True)
    X_test = test_df.drop(['SNP','CHR','L2','CHISQ'],axis=1)
    X_test = np.array(X_test.iloc[:,:])
    w_test = np.array(test_df['L2'])
    y_test = np.array(test_df['CHISQ']) 
    ypred = N*X_test.dot(tau)+1
   
    inv_w = 1/w_test 
    persnp_wsse = inv_w.dot((ypred-y_test)**2)/np.sum(inv_w)
    print('corrected pwsse',persnp_wsse)
    persnp = pd.DataFrame(data=[persnp_wsse],columns=['pwsse'])
    persnp.to_csv(outfile_prefix+'_pwsse',index=False,sep='\t')
    return

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--run_reg',action='store_true')
    parser.add_argument('--no_intersections',action='store_true')
    parser.add_argument('--methods',nargs='+',help='input list of methods for example: OLS Lasso Elnet Ridge')
    parser.add_argument('--outfile_prefix')
    parser.add_argument('--ld')
    parser.add_argument('--sumstat')
    parser.add_argument('--weight_ld')
    parser.add_argument('--leave_chrom',type=int)
    parser.add_argument('--compute_best_pwsse',action='store_true')
    parser.add_argument('--true_tau_file')
    parser.add_argument('--test_file')
    parser.add_argument('--varbeta_loss', action='store_true')
    parser.add_argument('--annot_prefix')
    parser.add_argument('--coef_file')
    parser.add_argument('--true_varbeta_file')
    parser.add_argument('--correct_pwsse',action='store_true')
    parser.add_argument('--rsid')
    parser.add_argument('--m550_rsid')
    parser.add_argument('--reg_weights',help='approx reciprocal of conditional variance function')
    parser.add_argument('--compute_varbeta',action='store_true')
    parser.add_argument('--alpha_file')
    parser.add_argument('--just_weights',action='store_true',help='break after the weights are stored')
    args = parser.parse_args()

    if args.run_reg:
        if args.no_intersections:
            run_reg(args)
    elif args.compute_best_pwsse:
        compute_best_pwsse(args) 
    elif args.varbeta_loss:
        compute_varbeta_loss(args)
    elif args.correct_pwsse:
        correct_pwsse(args)
    elif args.compute_varbeta:
        taudf = pd.read_csv(args.coef_file,delim_whitespace=True)
        tau = np.array(taudf.iloc[:,:])
        chrom = args.leave_chrom
        compute_varbeta(tau,chrom,args.annot_prefix,args.outfile_prefix,args.rsid,args.m550_rsid)
