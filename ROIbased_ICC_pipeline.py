import os
import pandas as pd
import ROIbased_ICC
import matplotlib.pyplot as plt


mzid = pd.read_csv('/nfs/t3/workingshop/liqinqin/projects/genetic_imaging/HCPdata/twin_data/random_50_mzid.csv',
                   index_col=0)
dzid = pd.read_csv('/nfs/t3/workingshop/liqinqin/projects/genetic_imaging/HCPdata/twin_data/random_50_dzid.csv',
                   index_col=0)
twinid = pd.concat((mzid, dzid))

stem_path = '/nfs/t3/workingshop/liqinqin/projects/genetic_imaging/'
vtx_mask_l = pd.read_csv(stem_path + 'data_analysis/searchlight/output/mask/roi/lh_all_mask.csv', index_col=0)
vtx_mask_r = pd.read_csv(stem_path + 'data_analysis/searchlight/output/mask/roi/rh_all_mask.csv', index_col=0)

vtx_map_l = list(
    pd.read_csv('/nfs/t3/workingshop/liqinqin/projects/genetic_imaging/data_analysis/mask/vertex_map_l.csv',
                index_col=0)['0'])
vtx_map_r = list(
    pd.read_csv('/nfs/t3/workingshop/liqinqin/projects/genetic_imaging/data_analysis/mask/vertex_map_r.csv',
                index_col=0)['0'])

colors_dic = {'body':['b', 'b', 'b', 'b', 'b'],
              'face':['g','g', 'g'],
              'place':['r','r', 'r', 'r'],
              'tool':['c','c', 'c', 'c']}


conditions = ['body', 'face', 'place', 'tool']
selective_dic = {'body': ['LOS', 'MTG', 'ITG', 'OTS'],
                 'face': ['OFA', 'FFA'],
                 'place': ['TOS', 'PPA', 'RSC'],
                 'tool': ['LO', 'ITGobject', 'pFs']}



all_regions = ['body','LOS', 'MTG', 'ITG', 'OTS',
                'face','OFA', 'FFA',
                'place','TOS', 'PPA', 'RSC',
                'tool','LO', 'ITGobject', 'pFs']

# gss_mask_stem_path = '/nfs/t3/workingshop/liqinqin/projects/genetic_imaging/data_analysis/finilize_icc/gss_mask/GSS_twin/thr0/'
# # thrs_folder = sorted([i for i in os.listdir(gss_mask_stem_path) if '.csv' not in i])
# avg_data_stem_path = '/nfs/t3/workingshop/liqinqin/projects/genetic_imaging/data_analysis/finilize_icc/'
# avg_output_path = os.path.join(avg_data_stem_path, 'regressout3con+gss+50twin+regresshm/')
# print(avg_output_path)

#1.calculate averaged_data

# input_path = '/nfs/t3/workingshop/liqinqin/projects/genetic_imaging/HCPdata/extractFromAWS/func/hp200_s4_level2_MSMALL.feat/beta/regressout3con/50twin/'
# ROIbased_ICC.gss_get_avg_data(input_path, gss_mask_stem_path, conditions, [vtx_map_l, vtx_map_r], twinid, avg_output_path )
# ROIbased_ICC.get_avg_data(input_path, [vtx_mask_l,vtx_mask_r], [vtx_map_l, vtx_map_r], all_regions, avg_output_path)
# ROIbased_ICC.get_topn_avg_data(input_path, [vtx_mask_l,vtx_mask_r], [vtx_map_l, vtx_map_r], all_regions, avg_output_path,npercent=0.4)

# rm outleirs
# ROIbased_ICC.rm_outliers_avg(avg_output_path, mzid, dzid)

#2.hist of averge data

# ROIbased_ICC.gss_avg_data_distribution(avg_output_path, selective_dic, colors_dic)
# ROIbased_ICC.gss_avg_data_distribution(avg_output_path+'rm_outliers/thr2/', selective_dic, colors_dic)

#3 sorted data into twin forms
# twin_form_path = ROIbased_ICC.sorted_avg_data_in_twin_forms(avg_output_path, mzid, dzid)
# twin_form_path = ROIbased_ICC.sorted_avg_data_in_twin_forms(avg_output_path+'rm_outliers/thr2/', mzid, dzid)
#
# # 4. sig icc masked resutls
#

icc_path = '/nfs/t3/workingshop/liqinqin/projects/genetic_imaging/data_analysis/finilize_icc/regressout3con+gss+50twin+regresshm/rm_outliers/thr2/results/'
ROIbased_ICC.sig_icc_results(icc_path)
#

# 5. descibe icc results

#
# ROIbased_ICC.des_icc(icc_results_path, conditions, selective_dic)
# ROIbased_ICC.des_cate_icc(icc_results_path, conditions, selective_dic)
#


#ploting
icc_path = '/nfs/t3/workingshop/liqinqin/projects/genetic_imaging/data_analysis/finilize_icc/regressout3con+gss+50twin+regresshm/rm_outliers/thr2/results/'
if not os.path.exists(icc_path+'plot/'):
    os.mkdir(icc_path+'plot/')

icc_files = sorted([i for i in os.listdir(icc_path) if '.csv' in i])


# icc_results = pd.DataFrame()
# for f in icc_files:
#     data = pd.read_csv(icc_path+f, index_col=0)
#     con = f.split('_')[0]
#     icc_results[con] = data[['iccmz','iccdz']].loc[con]
#
# icc_results.T.plot.bar(grid=True, width=0.65, color=['black','gray'])
# # plt.show()
# plt.savefig(icc_path+'plot/all_sys_mzdz_bar.png')

icc_results = pd.DataFrame()
regions = selective_dic['place']
# regions = ['body','face','place','tool']
for f in icc_files:
    data = pd.read_csv(icc_path+f, index_col=0)
    con = f.split('_')[0]
    icc_results[con] = data['mz-dz'].loc[regions]
# plt.show()
icc_results.plot(kind='bar', width=0.65, legend=False, grid=True)
plt.savefig(icc_path+'plot/cate_place.png')