# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 10:56:15 2017

@author: LQQ

@purpose: data preprosessing
"""

import numpy as np

import scipy.stats as stats
import pandas as pd
from ATT.algorithm import tools
from scipy.stats import ks_2samp
from scipy.stats import ttest_1samp

#%% data preparation 
class data_preparation(object):
    """
    this class is for data preparation including remove outliers, regress out covariates,and all analysis are based on pandas.Dataframe
    """
    def __init__(self,data):
        """

        :param data: data you want to preprocessed
        :return:
        """
        self.data = data.dropna() # drop nan data
        
    def regressout_signal(self,cov_sig):
        ori_data = np.array(self.data)[:,0:360]
        regre_data = np.array([tools.regressoutvariable(i,cov_sig) for i in ori_data.T]).T
        self.data = pd.DataFrame(regre_data,index = self.data.index,columns = self.data.columns[0:360]) 
        
        
    def all_roi_remove_outleirs(self,thrs = [-2,2]):
        """

        :param thrs: IQR threshold to remove outliers
        :return:num_of_outliers: cause all data is removed within the data frame, so only the number of removed outleirs will returen
        """
        
        num_of_outliers = []
        for i in self.data.columns:
            n,self.data[i] = tools.removeoutlier(self.data[i],meth = 'iqr',thr = thrs) #该方法在DataFrame内部修改值，所以不用返回值
            num_of_outliers.append(n)
        return num_of_outliers


    def twins_identifier(self,subjectID,twins_list):
        """
        :param:
            subjectID: an int or a list, containing the id of subjects whose twin-pair you want to find.
            twins_list: array-like matrix containing all the twins id. Its format must like this
    
                    twin1           twin2
                    -----------------------
                    10203           23234
                    20123           13454
                    23212           35433
                    ...
        :return:
            twin_pair: correspoding twin id of input twin list
        """
        twins_list = np.array(twins_list)
        twin1 = twins_list[:,0]
        twin2 = twins_list[:,1]
        
        if type(subjectID) == int:
            if subjectID in twin1:
                twin_pair = twin2[twin1==subjectID]
            elif subjectID in twin2:
                twin_pair = twin1[twin2 == subjectID]
            else:
                twin_pair = 'none'
    #        print('the twin pair of {} is {}'.format(subjectID,twin_pair))
        else:
            twin_pair = []
            i_pair = 0
            for i in subjectID:
                if i in twin1:
                    i_pair = twin2[twin1 == i]
                    twin_pair.append(int(i_pair))
                elif i in twin2:
                    i_pair = twin1[twin2 == i]
                    twin_pair.append(int(i_pair))
                else:
                    i_pair = 'none'
                    twin_pair.append(i_pair)
    #            print('the twin pair of {} is {}'.format(i,int(i_pair)))
        return twin_pair

    def rm_twin_pair(self,twinpair_id):
        """
        determine whehter a subject have its twin pair. If so, removing the twin pair and return the final number of removing outliers
        :param twinpair_id: a dataFrme contains informtaion of twin pair
        :return: num_of_final_outliers: number of total outliers when removing data of one subject's twin pair
        """

        num_of_final_outliers = []
        for i in self.data.columns:
            twin1_id = self.data[i][self.data[i].isnull()].index
            twin2_id = self.twins_identifier(twin1_id,twinpair_id)
            for j in twin2_id:
                if j != 'none':
                    self.data[i][j] = np.nan
            num_of_final_outliers.append(len(np.unique(list(twin1_id)+list(twin2_id)))) 
        return num_of_final_outliers    
    
    def save_data(self,output,filename):
        """
        :param output:output path
        :param filename: filename
        :return: none
        """
        self.data.to_csv(output+filename)
        
        

class twin_data_preparation(object):
    """
    1)This class in totally based on Dataframe;
    2)Purpose of this task is to organize data in twin's form
    """
    
    def __init__(self,twinid,in_twin_forms = True):
        """

        :param twinid: DataFrame contains twin infomation whose index is family ID
        :param in_twin_forms:Defaut True. If True, we organize twin data in two-columned form. Or we just extract twin data without any formation.
        :return:
        """
        self.twinid = twinid
        self.in_twin_forms = in_twin_forms
        self.data = pd.DataFrame()

        
        
    def load_univariate_twin_data(self,data,roi_name = 't', merge_to_twindata = False):
        """
        load univariate twin data
        :param data: data you want to load
        :param roi_name:  name of your variable, 't' by dufault
        :param merge_to_twindata: whether add to attribute(self.data) of this Class
        :return:
        """

        twinid1 = []
        twinid2 = []
        family_id = []
        for i in range(0,len(self.twinid)):
            if (self.twinid.ix[i][0] in np.array(data.index)) & (self.twinid.ix[i][1] in np.array(data.index)):
                twinid1.append(self.twinid.ix[i][0])
                twinid2.append(self.twinid.ix[i][1])
                family_id.append(self.twinid.index[i])
        
        if self.in_twin_forms == False:
            twin_data = data[data.index.isin(twinid1+twinid2)]
        elif self.in_twin_forms == True:
            if data.ndim == 1:
                twin1_data, twin2_data= data[twinid1],data[twinid2]
                twin_data = pd.DataFrame(np.array([twin1_data,twin2_data]).T,index = family_id,columns=[roi_name+'1',roi_name+'2'] )
            if data.ndim == 2:
                print(
                """we can not arrange data in twin forms when there is more than one column in the dataset
                Please change your data type to panda.series or set in_twin_form = Flase
                """)
                twin_data = None
        if merge_to_twindata == True:
            self.data[[roi_name+'1',roi_name+'2']] = twin_data
            self.screened_family_id = family_id
            return self.data
        else:
            return twin_data
    
    def load_covariates_for_twin(self,covdata,covname,merge_to_twindata = True,is_same = False):#this function require subject ID
        """
        add covariates once at a time to twin data
        parameters:
            covdata: dataframe contains two columns, one column should be the subject id (index in dataframe), and the 
                    other column is the value of that covairate
            covname: name of that covariates, used as column's name in Dataframe.
            merge_to_twindata: True or False. True means adding this covarite to self.data, while False means not.
            is_same: True or False. False by default. if you think value of covariate will be the same across twin
                    you can set this value to True. if not, this funcition will do comparision. 
                    if value of covariate are the same across twin order, one column of covariate will be added to the data,
                    if not, two columns will be added.
        """
        twinid1 = self.twinid.ix[self.data.index].iloc[:,0]
        twinid2 = self.twinid.ix[self.data.index].iloc[:,1]
        cov1 = covdata[covname].ix[twinid1]
        cov2 = covdata[covname].ix[twinid2]
        if (list(cov1.values) == list(cov2.values)) | (is_same == True):
            self.data[covname] = cov1.values
        else:
            self.data[covname+'1'] = cov1.values
            self.data[covname+'2'] = cov2.values#会按index插值吗？
        return self.data
        
    def load_sib_data(self,family_id_all,data,merge_to_twindata = True):
        
        family_id = list(self.data.index)
        
        sib_id = family_id_all[['s1','s2','s3','s4']].ix[family_id]
        
        sib_id['s1data'] = np.array(data.ix[sib_id['s1']])
        sib_id['s2data'] = np.array(data.ix[sib_id['s2']])
        sib_id['s3data'] = np.array(data.ix[sib_id['s3']])
        sib_id['s4data'] = np.array(data.ix[sib_id['s4']])

        sib_family = []
        sib_final_id = []
        sib_final_data = []
        
        for i in range(0,len(sib_id)):
            
            sib_id_data = sib_id.iloc[i,4:8]
            col = sib_id_data.dropna().index
            if len(col) > 0:
                sib_family.append(sib_id.iloc[i].name)
                sib_final_id.append(sib_id[col[0][0:2]].iloc[i])
                tar = col[0]
                sib_final_data.append(sib_id.iloc[i][tar]) 
                
        sib = pd.DataFrame(np.array([sib_final_id,sib_final_data]).T, index = sib_family,columns = [['sibID','sibData']])
        if merge_to_twindata == True:
            self.data = pd.concat([self.data,sib],axis = 1)
            self.sib_family = sib_family
            self.sibid = self.data['sibID']
            return self.data
        else:
            return sib
    
    def load_covariates_for_sib(self,covdata,covname,merge_to_twin_data = True): #should containing family data and 
        cov = covdata[covname].ix[self.sibid]
        if merge_to_twin_data == True:
            for i in range(0,len(covname)):
                self.data['sib_'+covname[i]] = np.array(cov[covname[i]])
            return self.data
        else:
            return cov 
        
    
    def save_data(self,output,filename):
        self.data.to_csv(output+filename)
        print('save '+filename+' successfully')
    
#%%
def isNormal(data):
    """
    whether a distribution is normal
    parameters:
        data: a dataFrame 
    return:
        a dataFrame contains test value and corresponding p value,and when p < 5%, distribution is normal
    """
    normaltest_results = []
    for i in data.columns:
        s,p = stats.normaltest(data[i].dropna())
        normaltest_results.append([s,p])
    return pd.DataFrame(normaltest_results,index = list(range(1,len(normaltest_results)+1)),columns=['statistic','p_value'])
 

def is_diff_distribution(data1,data2):
    """
    whether two distribution is the same
    parameters:
        data1: a dataFrame contains your data
        data2: a dataFrame contains your data, and it must have the same columns with data1
    return:
        diff_distribution: a dataFrame contains your statistic ans p value, and when p value > 5%, the two distribution is the same
    """
    diff_distribution = []
    for i in data1.columns:
        s,p = ks_2samp(data1[i],data2[i])
        diff_distribution.append([s,p])
    diff_distribution = pd.DataFrame(diff_distribution,index=data1.columns,columns=['statistic','pvalue'])
    return diff_distribution
    

def one_sample_ttest(data):
    """
    doing one sample t test to all rois 
    parameters:
        data: a dataFrame contains all your rois data
    return:
        a dataFrame contains t value and p value of each roi
    """
    roi_columns = data.columns
    results = []
    for i in roi_columns:
        t,p = ttest_1samp(data[i].dropna(),0)
        results.append([t,p])
    return pd.DataFrame(results,index = list(range(1,len(results)+1)),columns = ['t','p'])


def roi_correction(one_samp_test,thr):
    """
    correcting results of one-sample t test (reserving data have t > 0 )
    paremeters:
        one_samp_test: one_sample test resutls
        thr: threshold
    return:
        roi_id_pass_bonf: id of roi that pass bonf correction
        roi_id_pass_FDRbh: id of roi that pass FDRbh correction
    """
    pCorr = tools.PCorrection(np.array(one_samp_test)[0:360,1])
    roi_id_pass_FDRbh = one_samp_test.index[(one_samp_test.t>0)&(one_samp_test.p < pCorr.fdr_bh(thr))]
    roi_id_pass_bonf = one_samp_test.index[(one_samp_test.t>0)&(one_samp_test.p < pCorr.bonferroni(thr))]
    return roi_id_pass_bonf,roi_id_pass_FDRbh
    
def screen_data(input_data,ssid,vertex,save_out_file):
    screened_data = input_data.ix[ssid]
    screened_data = screened_data[vertex]
    return screened_data

