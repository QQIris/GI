import GI.funcs.ROIbased_pattern as ROIbased_pattern
import os
import pandas as pd
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


# 1.calculate pattern similirity
input_path = '/nfs/t3/workingshop/liqinqin/projects/genetic_imaging/HCPdata/extractFromAWS/func/hp200_s4_level2_MSMALL.feat/beta/regressout3con/50twin+regresshm/'
gss_mask_path = '/nfs/t3/workingshop/liqinqin/projects/genetic_imaging/data_analysis/finilize_icc/gss_mask/GSS_twin/thr0/'
pattern_simi_path = '/nfs/t3/workingshop/liqinqin/projects/genetic_imaging/data_analysis/finilize_pattern/nogss+50twin/'

# ROIbased_pattern.cal_pattern_simi_gss(input_path, gss_mask_path, conditions, [vtx_map_l, vtx_map_r], twinid, output_path)
# pattern_simi_path = ROIbased_pattern.cal_pattern_simi(input_path, [vtx_mask_l,vtx_mask_r], [vtx_map_l, vtx_map_r], all_regions, twinid, pattern_simi_path)

#2.display distribution of pattern similarity
# 3. describe heritability
# colors = ['b','g','r','c']
# input_path = '/nfs/t3/workingshop/liqinqin/projects/genetic_imaging/data_analysis/ROIbased4mm/original/pattern/similarity/gss_avg/'
# ROIbased_pattern.des_pattern_simi(input_path, selective_dic, colors)

# #3.sorted data for anova analyhsis
# ROIbased_pattern.sorted_data_for_anova_analysis(pattern_simi_path)
# #
# # input_path = '/nfs/t3/workingshop/liqinqin/projects/genetic_imaging/data_analysis/ROIbased4mm/original/pattern/similarity/unique/'
ROIbased_pattern.sorted_data_for_anova_plotting(pattern_simi_path)
#
# # ROIbased_pattern.sort_data_for_spss(pattern_simi_path)
#
# 4. three way anova plotting
input_path = pattern_simi_path + 'for_anova_plot/'
ROIbased_pattern.anova_plot(input_path, conditions, selective_dic)


#plot
# simi_path = '/nfs/t3/workingshop/liqinqin/projects/genetic_imaging/data_analysis/finilize_pattern/nogss+50twin+regresshm/'
# if not os.path.exists(simi_path+'des_plot/'):
#     os.mkdir(simi_path+'des_plot/')
# simi_files = sorted([i for i in os.listdir(simi_path) if '.csv' in i])
# simi_mean = pd.DataFrame()
# simi_sem = pd.DataFrame()
# for f in simi_files:
#     data = pd.read_csv(simi_path + f, index_col=0)
#     con = f.split('_')[0]
#     simi_mean[con] = data.groupby('twin')[con].mean()
#     simi_sem[con] = data.groupby('twin')[con].sem()
#
# simi_mean = simi_mean.sort_index(ascending=False)
# simi_sem = simi_sem.sort_index(ascending=False)
# simi_mean.T.plot.bar(grid=True, width=0.65, color=['black','gray'], yerr=simi_sem.T)
# # plt.show()
# plt.savefig(simi_path+'des_plot/sys.png')