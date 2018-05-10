import froi
from froi.algorithm.meshtool import get_n_ring_neighbor
# from froi.algorithm.tools import bfs
from scipy.stats.stats import pearsonr
import pyvttbl as pt
import nibabel as nib
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from scipy.stats import pearsonr

from ATT.algorithm import tools
import os
from nibabel.cifti2 import cifti2 as ci
import time

# def zTransform(r,n):
#     z = np.log((1+r)/(1-r)) * (np.sqrt(n-3)/2)
#     p = zprob(-z)
#     return p

def get_vertex_map(tdata_file):

    func_img = nib.load(tdata_file)
    h = func_img.header.get_index_map(1)
    tmp_l=[i for i in h.brain_models if i.brain_structure=='CIFTI_STRUCTURE_CORTEX_LEFT'][0]
    vertex_map_l = list(tmp_l.vertex_indices)
    tmp_r=[i for i in h.brain_models if i.brain_structure=='CIFTI_STRUCTURE_CORTEX_RIGHT'][0]
    vertex_map_r = list(tmp_r.vertex_indices)

    vertex_count = [tmp_l.index_count, tmp_r.index_count]

    return vertex_map_l,vertex_map_r,vertex_count


def pattern_similarity(data, mzid, dzid):

    n_vtx = len(data.columns)
    result = [list(pearsonr(data.ix[id1],data.ix[id2])) + [n_vtx,'mz'] for id1,id2 in mzid.values] # pattern similarity for mz twins
    result +=  [list(pearsonr(data.ix[id1],data.ix[id2])) + [n_vtx,'dz'] for id1,id2 in dzid.values] # pattern similarity for dz twins

    result = pd.DataFrame(result,index=list(mzid.index)+list(dzid.index),columns=['r','p','num_vtx','twin_type'])
    result['fisher_z'] = tools.r2z(np.array(result['r']))

    mean_mz = np.mean(result['fisher_z'][result.twin_type == 'mz']) # mean of pattern similarity
    mean_dz = np.mean(result['fisher_z'][result.twin_type == 'dz']) # mean of pattern similarity
    ttest = ttest_ind(result['fisher_z'][result.twin_type == 'mz'],result['fisher_z'][result.twin_type == 'dz'])

    num_of_valid_mz = sum([1 for i in result['p'][result.twin_type == 'mz'] if i > 0.05])

    diff_result = [mean_mz,mean_dz,mean_mz-mean_dz,len(data.columns),ttest[0],ttest[1],num_of_valid_mz]
    return diff_result,list(result['fisher_z'])

def level_analysis_vtx(data, masks_of_interest, mzid, dzid,ignore_hemi=False):
    """

    :param data: a list if ignore_hemi==False, contain data of left hemi and right hemi;
    :param masks_of_interest: a list if ignore hemi == False ,contiaon mask for leftcan be n-ring masks of given vtx or masks of ROI.index of two masks must be the same.
    :param mzid:
    :param dzid:
    :return:
    """
    level_result = []
    detail_result =[]

    if ignore_hemi == False:
        for i in list(masks_of_interest.index):
            print(i)
            mask = masks_of_interest.ix[i] # get mask for each vertex
            mask = mask[~mask.isnull()] # remove NaN
            diff_result,result = pattern_similarity(data[mask], mzid, dzid)
            level_result.append(diff_result)
            detail_result.append(result)

        level_result = pd.DataFrame(level_result, index=masks_of_interest.index, columns=['meanMZ', 'meanDZ', 'meanMZ-meanDZ', 'vtx_num', 't', 'p', 'num_of_valid_mz_005']) # for level analysis
        detail_result = pd.DataFrame(detail_result, index=masks_of_interest.index, columns=list(mzid.index) + list(dzid.index)).T # for ANOVA analysis

    if ignore_hemi == True:
        data_l = data[0]
        data_r = data[1]
        mask_for_l = masks_of_interest[0]
        mask_for_r = masks_of_interest[1]
        for i in list(mask_for_l.index):
            print(i)
            mask_l = mask_for_l.ix[i]
            mask_r = mask_for_r.ix[i]
            mask_l = mask_l[~mask_l.isnull()]
            mask_r = mask_r[~mask_r.isnull()]
            data_used = pd.concat((data_l[mask_l],data_r[mask_r]),axis=1)
            diff_result,result = pattern_similarity(data_used,mzid,dzid)
            level_result.append(diff_result)
            detail_result.append(result)

        level_result = pd.DataFrame(level_result, index=mask_for_l.index, columns=['meanMZ', 'meanDZ', 'meanMZ-meanDZ', 'vtx_num', 't', 'p', 'num_of_valid_mz_005']) # for level analysis
        detail_result = pd.DataFrame(detail_result, index=mask_for_l.index, columns=list(mzid.index) + list(dzid.index)).T # for ANOVA analysis

    return level_result,detail_result



def anova_analysis(data,dv,wfactors,bfactors,alpha=0.05,return_type='all'):
    """

    :param data:
    :param dv: dependent variables
    :param wfactors: list of within subject variables
    :param bfactors: list of between subject variables
    :param alpha:
    :param return_type: the results you concerned about. 'main' for main effect, 'interaction' for interaction effect, 'all' for all.
    :return: A dataFrame contain 'F','ci','mse','eta','p' for each vertex
    """
    df = pt.DataFrame()
    df['subject_id'] = data.index
    df[wfactors[0]] = data[wfactors[0]]
    df[bfactors[0]] = data[bfactors[0]]
    aov_result = []
    values = ['F','ci','mse','eta','p']
    for column in data.columns[:-2]:
        print('anova for ',column)
        df[dv] = data[column]
        aov = df.anova(dv,sub='subject_id',wfactors=wfactors,bfactors=bfactors,alpha=alpha)

        if return_type == 'main':
            keys = aov.keys()[0:2]

        elif return_type == 'interaction':
            keys = aov.keys()[2]
        elif return_type == 'all':
            pass

        aov_result.append([aov[keys][value] for value in values])

    aov_result = pd.DataFrame(aov_result,index=data.columns[:-2],columns=values)

    return aov_result

def to_cifti_file(data_file,vertex_map,write_img,output_name,columns='p'):
    """

    :param data_file:a list of data path, each element in it must be a list contain left hemi data file and right hemi data path
    :param vertex_map:
    :param write_img:
    :param output_name:
    :param columns:
    :return:
    """


    data_l = pd.read_csv(data_file[0],index_col=0)[columns].to_dict()
    data_r = pd.read_csv(data_file[1],index_col=0)[columns].to_dict()

    vtx_l = data_l.keys()#get vertex id of data
    vtx_r = data_r.keys()

    write_data = np.zeros_like(write_img.get_data())
    for vtx in vtx_l:
        write_data[0][0:29696][np.array(vertex_map[0]) == vtx] = data_l[vtx] #write data according to vertex_map
    for vtx in vtx_r:
        write_data[0][29696:59412][np.array(vertex_map[1]) == vtx] = data_r[vtx]

    new_img = ci.Cifti2Image(write_data,write_img.header,write_img.nifti_header)
    ci.save(new_img,output_name)
    print('successfully save CIFTI file for data ', data_file)

def level_specific_searchlight_pipeline(data_file,correspond_condition,vtx_mask,tfile,output_path,mzid,dzid):
    """
    based on the given file do level heritability analysis
    :param data_file: a string or a list of string that contains path of the data file (csv fileformat)
    :param correspond_condition: list of correspond condition name of your data file
    :param vtx_mask:
    :param tfile: tfile to get vtx_map of each column of input data
    :param output_path: path to save results
    :return:
    """
    begin = time.time()

    vtx_map_l,vtx_map_r,vtx_count = get_vertex_map(tfile)
    mask_l = vtx_mask[0]
    mask_r = vtx_mask[1]

    for i,file in enumerate(data_file):
        data = pd.read_csv(file,index_col=0)
        print('-----------------------------------------------------------------------')
        print(file)
        data.columns = [vtx_map_l+vtx_map_r] #set columns as vertex id
        data_l = data.iloc[:,:vtx_count[0]]#seperate data into right hemi and left hemi
        data_r = data.iloc[:,vtx_count[0]:]

        level_result_r,detail_result_r = level_analysis_vtx(data_r,mask_r,mzid,dzid)
        level_result_l,detail_result_l = level_analysis_vtx(data_l,mask_l,mzid,dzid)

        level_result_l.to_csv(output_path+'lh_'+correspond_condition[i]+'searchlight_MVPA_similarity_diff_of_MZDZ_mean_result.csv')
        level_result_r.to_csv(output_path+'rh_'+correspond_condition[i]+'searchlight_MVPA_similarity_diff_of_MZDZ_mean_result.csv')
        detail_result_l.to_csv(output_path+'lh_'+correspond_condition[i]+'searchlight_MVPA_similarity_diff_of_MZDZ_fisherz_for_each_subject.csv')
        detail_result_r.to_csv(output_path+'rh_'+correspond_condition[i]+'searchlight_MVPA_similarity_diff_of_MZDZ_fisherz_for_each_subject.csv')
        print('running time for one data file: ', time.time()-begin)


def criteria_for_category(data_file,condition_file,non_condition_file):
    data = pd.read_csv(data_file,index_col=0)
    condition = pd.read_csv(condition_file,index_col=0)[['p','meanMZ-meanDZ']]
    non_condition = pd.read_csv(non_condition_file,index_col=0)
    screened_mask = [1 if i<0.05 else 0 for i in condition['p'].values] # criteria1: heritability of condition must be significant
    np.array(screened_mask)[condition['meanMZ-meanDZ'] <= non_condition['meanMZ-meanDZ']] = 0 #criteria2: effective size of one condition must larger than non-condition

    data['screened_F'] = data['F'] * screened_mask
    data['screened_p'] = data['p'] * screened_mask

    data.to_csv(data_file)
    print('suceesfully processed ', data_file.split('/')[-1])


def Intraclass_correlation(data,mzid,dzid):
    icc_results = []
    print('now we are performing ICC analysis')

    for i in list(data.columns):
        print(i)
        mzdata1 = data[i].ix[mzid['twin1']]
        mzdata2 = data[i].ix[mzid['twin2']]
        dzdata1 = data[i].ix[dzid['twin1']]
        dzdata2 = data[i].ix[dzid['twin2']]
        rmz,pmz = pearsonr(mzdata1.dropna(),mzdata2.dropna())
        rdz,pdz = pearsonr(dzdata1.dropna(),dzdata2.dropna())
        icc_results.append([rmz,pmz,rdz,pdz,rmz-rdz])

    icc_results = pd.DataFrame(icc_results,index=[i[0] for i in data.columns],columns=['rmz','pmz','rdz','pdz','rmz-rdz'])
    mask_mz_005 = [1 if i<0.05 else 0 for i in icc_results['pmz'].values]
    mask_dz_005 = [1 if i<0.05 else 0 for i in icc_results['pdz'].values]
    icc_results['sig005mz'] = icc_results['rmz'] * mask_mz_005
    icc_results['sig005dz'] = icc_results['rdz'] * mask_dz_005
    icc_results['sigmz-dz'] = icc_results['rmz-rdz'] * mask_mz_005
    return  icc_results


def create_final_mask(hemi,anova_path,level_path,anova_files,level_files,conditions,output,thr=0.05):


    print(hemi)
    final_mask = pd.DataFrame()

    for condition in conditions:
        level_data = pd.read_csv(level_path+[i for i in level_files if (condition in i) & (hemi in i)][0],index_col=0)
        level_compare_data_files = [i for i in level_files if (hemi in i) & (condition not in i) & ('avgall' not in i)]
        c_data1 =  pd.read_csv(level_path + level_compare_data_files[0],index_col=0)['meanMZ-meanDZ']
        c_data2 =  pd.read_csv(level_path + level_compare_data_files[1],index_col=0)['meanMZ-meanDZ']
        c_data3 =  pd.read_csv(level_path + level_compare_data_files[2],index_col=0)['meanMZ-meanDZ']
        anova_con_files = [i for i in anova_files if (condition in i) & (hemi in i)]
        a_data1 = pd.read_csv(anova_path+anova_con_files[0],index_col=0)['p']
        a_data2 = pd.read_csv(anova_path+anova_con_files[1],index_col=0)['p']
        a_data3 = pd.read_csv(anova_path+anova_con_files[2],index_col=0)['p']


        mask1,mask2,mask3,mask4 = [np.zeros_like(c_data1),np.zeros_like(c_data1),np.zeros_like(c_data1),np.zeros_like(c_data1)]
        mask1[np.array(level_data['meanMZ-meanDZ'] > c_data1) & np.array(level_data['meanMZ-meanDZ'] > c_data2) & np.array(level_data['meanMZ-meanDZ'] > c_data3)] = 1 #signifant effective size than other region
        mask2[np.array(a_data1 < thr) & np.array(a_data2 < thr) & np.array(a_data3 < thr)] = 1 #significant interaction effect in each conditon
        # mask3[np.array(anova_non<0.01)] = 1 # significant condtion vs non-condition effect
        mask4[np.array(level_data['p'] < thr)] = 1 #significant heritability within an conditon
        mask = mask1*mask2*mask4
        print(condition)
        print(sum(mask))

        final_mask[condition[1:]] = mask

    final_mask.index = a_data1.index
    final_mask.to_csv(anova_path+'final_mask/'+hemi+output)

def prepare_anova_data():
    pass