import GI.searchlight
import pandas as pd
import time
import os
import nibabel as nib
import numpy as np

stem_path = '/nfs/t3/workingshop/liqinqin/projects/genetic_imaging/'
mask_path = 'data_analysis/domain_specific_roibased/mask/temporal_cortex/'
input_path = 'data_analysis/searchlight_level_vs_category/input/beta/'
output_path_level = 'data_analysis/searchlight_level_vs_category/output/beta_results/level/'
output_path_cate = 'data_analysis/searchlight_level_vs_category/output/beta_results/category/'
tfile = '/nfs/t3/workingshop/liqinqin/projects/genetic_imaging/HCPdata/others/tstat1.dtseries.nii'


mzid= pd.read_csv('/nfs/t3/workingshop/liqinqin/projects/genetic_imaging/HCPdata/twin_data/mzid_s8.csv',index_col=0)
dzid = pd.read_csv('/nfs/t3/workingshop/liqinqin/projects/genetic_imaging/HCPdata/twin_data/dzid_s8.csv',index_col=0)


vtx_map_l,vtx_map_r,vtx_count = searchlight.get_vertex_map(tfile)
vtx_mask_l = pd.read_csv('/nfs/t3/workingshop/liqinqin/projects/genetic_imaging/data_analysis/searchlight_level_vs_category/output/mask/vtx/lh_whole_4ring_neighbor_mask_cleaned.csv',index_col = 0)
vtx_mask_r = pd.read_csv('/nfs/t3/workingshop/liqinqin/projects/genetic_imaging/data_analysis/searchlight_level_vs_category/output/mask/vtx/rh_whole_4ring_neighbor_mask_cleaned.csv',index_col = 0)


# ----------------------------------------------------level_analysis---------------------------------------------------------------------------
#
# begin = time.time()
#
# mask_l = vtx_mask_l
# mask_r = vtx_mask_r
# data_stem = '/nfs/t3/workingshop/liqinqin/projects/genetic_imaging/HCPdata/extractFromAWS/func/hp200_s12_level2.feat/'
# data_file = sorted([i for i in os.listdir(data_stem) if ('fix' not in i) & ('_avg' not in i)])



# file = data_file[0]
# for file in data_file:
#
# condition = file[:-26]
#
# data = pd.read_csv(data_stem+file,index_col=0)
# print('-----------------------------------------------------------------------')
# print(file)
# print(condition)
#
# data.columns = [vtx_map_l+vtx_map_r] #set columns as vertex id
# data_l = data.iloc[:,:vtx_count[0]]#seperate data into right hemi and left hemi
# data_r = data.iloc[:,vtx_count[0]:]


# #not merge hemi
# output_path = stem_path+output_path_level+'4ring/s12/'
# # output_path = stem_path + 'data_analysis/searchlight_level_vs_category/test_with_random_twin/pattern/s4/'
# level_result_r,detail_result_r = searchlight.level_analysis_vtx(data_r,mask_r,mzid,dzid)
# level_result_l,detail_result_l = searchlight.level_analysis_vtx(data_l,mask_l,mzid,dzid)
# level_result_l.to_csv(output_path+'lh_'+condition+'_4ring_similarity_diff_of_MZDZ_mean_result.csv')
# level_result_r.to_csv(output_path+'rh_'+condition+'_4ring_similarity_diff_of_MZDZ_mean_result.csv')
# detail_result_l.to_csv(output_path+'detail/lh_'+condition+'_4ring_similarity_diff_of_MZDZ_fisherz_for_each_subject.csv')
# detail_result_r.to_csv(output_path+'detail/rh_'+condition+'_4ring_similarity_diff_of_MZDZ_fisherz_for_each_subject.csv')
#
# print(time.time()-begin)

#
# #merge hemi
# #
# # level_result,detail_result = searchlight.level_analysis_vtx([data_l,data_r],[vtx_mask_l,vtx_mask_r],mzid,dzid,ignore_hemi=True)
# # level_result.to_csv(output_path+condition+'_pvc_merge_hemi_similarity_diff_of_MZDZ_mean_result.csv')
# # detail_result.to_csv(output_path+'detail/'+condition+'_pvc_merge_hemi_similarity_diff_of_MZDZ_mean_result.csv')
# # print(time.time()-begin)
#
# # ICC analysis

# output_path = stem_path+output_path_level+'ICC_vertex/s12/'
# # output_path = stem_path + 'data_analysis/searchlight_level_vs_category/test_with_random_twin/ICC/s4/'
# icc_l = searchlight.Intraclass_correlation(data_l,mzid,dzid)
# icc_l.to_csv(output_path+'lh_'+condition+'_diff_of_icc_mzdz_results_of_each_vertex.csv')
# icc_r = searchlight.Intraclass_correlation(data_r,mzid,dzid)
# icc_r.to_csv(output_path+'rh_'+condition+'_diff_of_icc_mzdz_results_of_each_vertex.csv')
# print(time.time()-begin)


# # #-----------------------------------anova analysis-----------------------------------------------------------------------------------

# prepare data
#
# level_path = '/nfs/t3/workingshop/liqinqin/projects/genetic_imaging/data_analysis/searchlight_level_vs_category/output/beta_results/level/4ring/s12/detail/'
# anova_path = stem_path + input_path + '4ring/s12/for_anova/unique/'
# condition = ['_body','_face','_place','_tool']
# # non_condition = ['_nonbody','_nonface','_nonplace','_nontool']
# twin = ['MZ' for i in range(0,71)]
# twin += ['DZ' for i in range(0,41)]
#
#
# files = sorted([i for i in os.listdir(level_path) for j in condition if j in i])
# # com_files  = sorted([i for i in os.listdir(level_path) for j in non_condition if j in i])
# # file_set = [[files[i],com_files[i]] for i in range(0,8)]
#
# file_set = [[files[i],files[j]] for i in range(0,4) for j in range(0,4) if i < j]
# file_set += [[files[i],files[j]] for i in range(4,8) for j in range(4,8) if i < j]
#
# for i in file_set:
#     hemi = i[0].split('_')[0]
#     condition = i[0].split('_')[1]
#     print(condition)
#     com_comdition = i[1].split('_')[1]
#     data = pd.read_csv(level_path+i[0],index_col=0)
#     comdata = pd.read_csv(level_path+i[1], index_col=0)
#     data['twin'] = twin
#     comdata['twin'] = twin
#     data['condition'] = [condition for index in range(0,112)]
#     comdata['condition'] = [com_comdition for index in range(0, 112)]
#     final = pd.concat((data, comdata))
#     final.to_csv(anova_path+hemi+'_'+condition+'_vs_'+com_comdition+'_unique_beta_s12_fisherz_for_anova.csv')

# anova analysis
#
# begin = time.time()
# anova_path = stem_path + input_path + '4ring/s12/for_anova/unique/'
# data_file = sorted([i for i in os.listdir(anova_path)])[11]
# # for data_file in sorted([i for i in os.listdir(anova_path)]):
#
# print(data_file)
# data = pd.read_csv(anova_path+data_file,index_col=0)
# inter_result = searchlight.anova_analysis(data,'fisherz',['condition'],['twin'],return_type='interaction')
# print(inter_result)
# inter_result.to_csv(stem_path+output_path_cate+'4ring/s12/'+data_file[:-22]+'_condition_twin_interaction_anova_results.csv')
#
# print(time.time()-begin)

# pvc vs vvc
# region = 'body'
#
# data_file = stem_path + input_path + 'system_pvc_anova /' + region+'_unique_system_pvc_twin_for_anova.csv'
# data = pd.read_csv(data_file,index_col=0)
#
# rois = ['v1','v2','v3','v4']
# for i in rois:
#     tmp = data[(data['roi'] == i)|(data['roi'] == region)]
#     aov_result = searchlight.anova_analysis(tmp,'fisher_z',['roi'],['twin'],return_type='interaction')
#     print(aov_result)


# # #------------------------------------------------save to cifti file----------------------------------------------------------------------
begin = time.time()
write_img = nib.load(tfile)
path =  '/nfs/t3/workingshop/liqinqin/projects/genetic_imaging/data_analysis/searchlight_level_vs_category/output/beta_results/level/ICC_vertex/s8/'

files = [path +i for i in sorted(os.listdir(path)) if 'avg' in i]
file_set = [[files[0],files[5]],[files[1],files[6]],[files[2],files[7]],[files[3],files[8]],[files[4],files[9]]]
# file_set = [[path+'lh_final_mask_of_4ring_unique_of_anova_results_thre_001.csv',path+'rh_final_mask_of_4ring_unique_of_anova_results_thre_001.csv']]
# file_set = [[files[0],files[1]]]

for data_file in file_set:
    print('--------------------------------------------')
    print('save cifti files')
    print(data_file[0].split('/')[-1])

#     # icc resutls
    searchlight.to_cifti_file(data_file,[vtx_map_l,vtx_map_r],write_img,stem_path+'ICC_vertex/s8/cifti/mz/'+data_file[0].split('/')[-1][3:-44]+'_s4_msmall_sig_rmz_icc_value.dtseries.nii',columns='sig005mz')
    searchlight.to_cifti_file(data_file,[vtx_map_l,vtx_map_r],write_img,stem_path+'ICC_vertex/s8/cifti/dz/'+data_file[0].split('/')[-1][3:-44]+'_s4_msmall_sig_rdz_icc_value.dtseries.nii',columns='sig005dz')
    searchlight.to_cifti_file(data_file,[vtx_map_l,vtx_map_r],write_img,stem_path+'ICC_vertex/s8/cifti/mz_dz/'+data_file[0].split('/')[-1][3:-44]+'_s4_msmall_sigMZ_DZ_icc_value.dtseries.nii',columns='sigmz-dz')
#
# #     #general level resutls
#     searchlight.to_cifti_file(data_file,[vtx_map_l,vtx_map_r],write_img,stem_path+'data_analysis/searchlight_level_vs_category/test_with_random_twin/pattern/s2/cifti/'+data_file[0].split('/')[-1][3:-46]+'_beta_s2_msmall_sig_effective_size_value.dtseries.nii',columns='meanMZ-meanDZ')
#     searchlilght.to_cifti_file(data_file,[vtx_map_l,vtx_map_r],write_img,stem_path+'data_analysis/searchlight_level_vs_category/test_with_random_twin/pattern/s2/cifti/'+data_file[0].split('/')[-1][3:-46]+'_beta_s2_msmall_sig_p_value.dtseries.nii',columns='p')
# #
# #     #anova resutls
# #     # searchlight.to_cifti_file(data_file,[vtx_map_l,vtx_map_r],write_img,stem_path+output_path_cate+'4ring/s4/cifti/'+data_file[0].split('/')[-1][3:-44]+'_beta_s4_msmall_F_value.dtseries.nii',columns='F')
# #     # searchlight.to_cifti_file(data_file,[vtx_map_l,vtx_map_r],write_img,stem_path+output_path_cate+'4ring/s4/cifti/'+data_file[0].split('/')[-1][3:-44]+'_beta_s4_msmall_p_value.dtseries.nii',columns='p')
# #
#     searchlight.to_cifti_file(data_file,[vtx_map_l,vtx_map_r],write_img,stem_path+output_path_cate+'4ring/s12/final_mask/'+data_file[0].split('/')[-1][3:-45]+'4ring_for_body_thr005.dtseries.nii',columns='body')
#     searchlight.to_cifti_file(data_file,[vtx_map_l,vtx_map_r],write_img,stem_path+output_path_cate+'4ring/s12/final_mask/'+data_file[0].split('/')[-1][3:-45]+'4ring_for_face_thr005.dtseries.nii',columns='face')
#     searchlight.to_cifti_file(data_file,[vtx_map_l,vtx_map_r],write_img,stem_path+output_path_cate+'4ring/s12/final_mask/'+data_file[0].split('/')[-1][3:-45]+'4ring_for_place_thr005.dtseries.nii',columns='place')
#     searchlight.to_cifti_file(data_file,[vtx_map_l,vtx_map_r],write_img,stem_path+output_path_cate+'4ring/s12/final_mask/'+data_file[0].split('/')[-1][3:-45]+'4ring_for_tool_thr005.dtseries.nii',columns='tool')
#
# print(time.time()-begin)

#
# # #------------------------------------------------create final mask----------------------------------------------------------------------
#
# anova_path = stem_path + output_path_cate + '4ring/s12/'
# level_path = stem_path + output_path_level + '4ring/s12/'
# anova_files = sorted([i for i in os.listdir(anova_path) if 'unique' in i])
# level_files = sorted([i for i in os.listdir(level_path) if 'unique' in i])
# # level_data = None
# conditions = ['_body','_face','_place','_tool']
#
#
# output = 'beta_s12_final_mask_of_4ring_unique_of_anova_results_thre_005.csv'
# searchlight.create_final_mask('lh_',anova_path,level_path,anova_files,level_files,conditions,output)
# searchlight.create_final_mask('rh_',anova_path,level_path,anova_files,level_files,conditions,output)
#


#----------------------------------------define a pipeline--------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
#get nring neighbor
    # vtc_l_vtx = nib.freesurfer.read_label(stem_path+mask_path+'vtc_mask_l.label')
    # vtc_r_vtx = nib.freesurfer.read_label(stem_path+mask_path+'vtc_mask_r.label')
    #
    # img_l = nib.load('/nfs/t3/workingshop/liqinqin/projects/genetic_imaging/HCPdata/S3download/100307.L.white_MSMAll.32k_fs_LR.surf.gii')
    # vtc_l_mask = np.zeros(len(limg_l.darrays[0].data))
    # vtc_l_mask[vertex_map_l] = 1
    #
    # img_r = nib.load('/nfs/t3/workingshop/liqinqin/projects/genetic_imaging/HCPdata/S3download/100307.R.white_MSMAll.32k_fs_LR.surf.gii')
    # vtc_r_mask = np.zeros(len(img_r.darrays[0].data))
    # vtc_r_mask[vertex_map_r] = 1
    #
    #
    # vtc_3ring_l = get_n_ring_neighbor(img_l.darrays[1].data,3,mask = vtc_l_mask)
    # vtc_3ring_r = get_n_ring_neighbor(img_r.darrays[1].data,3,mask = vtc_r_mask)
    #
    # pd.DataFrame([vtc_3ring_l[i] for i in vtc_l_vtx],index=vtc_l_vtx).to_csv(stem_path+input_path+'vtc/vtc_3ring_l_mask.csv')
    # pd.DataFrame([vtc_3ring_r[i] for i in vtc_r_vtx],index=vtc_r_vtx).to_csv(stem_path+input_path+'vtc/vtc_3ring_r_mask.csv')



#-----------------------------ploting----------------------------------------------------------------------------------
# import matplotlib.pyplot as plt
#
# lh_vtx = pd.read_csv('/nfs/t3/workingshop/liqinqin/projects/genetic_imaging/data_analysis/searchlight_level_vs_category/output/mask/lh_all_mask_conbined.csv',index_col=0)
# rh_vtx = pd.read_csv('/nfs/t3/workingshop/liqinqin/projects/genetic_imaging/data_analysis/searchlight_level_vs_category/output/mask/rh_all_mask_conbined.csv',index_col=0)
#
#
# plot_path = '/nfs/t3/workingshop/liqinqin/projects/genetic_imaging/data_analysis/searchlight_level_vs_category/output/beta_results/level/'
# data_path = '4ring/s8/'
# files = sorted([plot_path+data_path+i for i in os.listdir(plot_path+data_path) if 'surface' in i])
#
# file_set = [[files[0],files[5]],[files[1],files[6]],[files[2],files[7]],[files[3],files[8]],[files[4],files[9]]]
#
#
# for f in file_set:
#     conditon = f[0].split('/')[-1][3:-44]
#     condition_name = conditon[:-25]
#
#     if condition_name == 'avg_four_co':
#         condition_name = 'system'
#
#     roi = ['v1','v2','v3','v4',condition_name]
#     data_l = pd.read_csv(f[0],index_col=0)
#     data_r = pd.read_csv(f[1],index_col=0)
#
#     lh_data = []
#     rh_data = []
#     whole_data = []
#     for i in lh_vtx.index:
#         lh_data.append(list(data_l['meanMZ-meanDZ'].ix[lh_vtx.ix[i].dropna()]))
#         rh_data.append(list(data_r['meanMZ-meanDZ'].ix[rh_vtx.ix[i].dropna()]))
#         whole_data.append(list(data_r['meanMZ-meanDZ'].ix[rh_vtx.ix[i].dropna()])+list(data_l['meanMZ-meanDZ'].ix[lh_vtx.ix[i].dropna()]))
#
#
#     lh_data = pd.DataFrame(lh_data,index=lh_vtx.index)
#     rh_data = pd.DataFrame(rh_data,index=lh_vtx.index)
#     whole_data = pd.DataFrame(whole_data,index=lh_vtx.index) # to pandas DataFrame
#
#
#     lh_data.ix[roi].mean(axis=1).plot(kind='bar',yerr = lh_data.ix[roi].sem(axis=1),fontsize=20,grid=True)
#     plt.savefig(plot_path+data_path+'plot/lh_'+condition_name+'_all_rois_mean_value_of_effective_size_on_pattern_level_analysis.png')
#     plt.close()
#
#     rh_data.ix[roi].mean(axis=1).plot(kind='bar',yerr = rh_data.ix[roi].sem(axis=1),fontsize=20,grid=True)
#     plt.savefig(plot_path+data_path+'plot/rh_'+condition_name+'_all_rois_mean_value_of_effective_size_on_pattern_level_analysis.png')
#     plt.close()
#
#     whole_data.ix[roi].mean(axis=1).plot(kind='bar',yerr = whole_data.ix[roi].sem(axis=1),fontsize=20,grid=True)
#     plt.savefig(plot_path+data_path+'plot/whole_'+condition_name+'_all_rois_mean_value_of_effective_size_on_pattern_level_analysis.png')
#     plt.close()