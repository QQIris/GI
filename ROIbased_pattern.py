import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats.stats import pearsonr


def cal_pattern_simi_gss(input_path, gss_mask_input, conditions, vtx_map, twinid, output_path):
    print('cal pattern similarity under gss mask')
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    beta_files = sorted([i for i in os.listdir(input_path) if '.csv' in i])
    gss_files = sorted([i for i in os.listdir(gss_mask_input) if '.csv' in i])
    gss_file_set = [[gss_files[i], gss_files[i + len(gss_files) / 2]] for i in range(0, len(gss_files) / 2)]
    rois = [i.split('_')[1] for i in gss_files[:len(gss_files) / 2]]
    twin = ['mz' if i < 50 else 'dz' for i in range(0, 100)]
    for con in conditions:
        f = [i for i in beta_files if con in i][0]
        print('-------------------------------------')
        print(f)
        data = pd.read_csv(input_path + f, index_col=0)
        data_l = data.iloc[:, :29696]
        data_r = data.iloc[:, 29696:59412]
        data_l.columns = vtx_map[0]
        data_r.columns = vtx_map[1]
        pattern_simi_detail = pd.DataFrame(index=twinid.index, columns=rois)
        for sets in gss_file_set:
            pattern_simi_temp = pd.DataFrame(index=twinid.index, columns=rois)
            gss_mask_l = pd.read_csv(gss_mask_input + sets[0], index_col=0)
            gss_mask_r = pd.read_csv(gss_mask_input + sets[1], index_col=0)
            roi = sets[0].split('_')[1]
            print(roi)
            for family in twinid.index:
                pattern_l = data_l.ix[list(twinid.ix[family])][gss_mask_l.ix[family].dropna()]
                pattern_r = data_r.ix[list(twinid.ix[family])][gss_mask_r.ix[family].dropna()]
                pattern = pd.concat((pattern_l, pattern_r), axis=1)

                pattern_simi_detail[roi].ix[family] = pattern.T.corr().values[0][1]
        pattern_simi_detail['twin'] = twin
        pattern_simi_detail['condition'] = [con for i in range(0, 100)]
        pattern_simi_detail.to_csv(output_path + '_'.join(f.split('_')[:-1]) + '_pattern_similirity_detail.csv')
        print('successful save csv file')


def cal_pattern_simi(input_path, vtx_mask, vtx_map, all_regions, twinid, output_path):
    print('calculate pattern simialrity without gss mask')
    files = sorted([i for i in os.listdir(input_path) if '_fix_' in i])
    pattern_simi_detail = pd.DataFrame(index=twinid.index)
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    for f in files:
        print(f)
        con = f.split('_')[0]
        data = pd.read_csv(input_path + f, index_col=0)
        data_l = data.iloc[:, :29696]
        data_r = data.iloc[:, 29696:59412]
        data_l.columns = vtx_map[0]
        data_r.columns = vtx_map[1]

        for region in all_regions:
            print(region)
            pattern_l = data_l[vtx_mask[0].loc[region].dropna()]
            pattern_r = data_r[vtx_mask[1].loc[region].dropna()]
            pattern = pd.concat((pattern_l, pattern_r), axis=1)
            pattern_simi_detail[region] = [pearsonr(pattern.loc[twin1], pattern.loc[twin2])[0] for twin1, twin2 in
                                           twinid.values]

        pattern_simi_detail['twin'] = ['mz' if i < 50 else 'dz' for i in range(0, len(twinid))]
        pattern_simi_detail['condition'] = [con for i in range(0, len(twinid))]

        pattern_simi_detail.to_csv(
            os.path.join(output_path, '_'.join(f.split('_')[:-1]) + '_pattern_similarity_detail.csv'))
    return output_path


def dstbs_of_pattern_simi(input, output):
    pass


def des_pattern_simi(input_path, selective_dic, colors):
    files = sorted([i for i in os.listdir(input_path) if '.csv' in i])
    output_path = os.path.join(input_path, 'describe')
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    for n, f in enumerate(files):
        con = f.split('_')[0]
        print(con)
        regions = [con] + selective_dic[con]
        data = pd.read_csv(input_path + f, index_col=0)
        data.groupby('twin')[regions].mean().plot(kind='bar', subplots=True, layout=(1, len(regions)),
                                                  sharey=True, legend=False, color=colors[n], width=1, ylim=(0, 0.5))
        plt.subplots_adjust(hspace=0, wspace=0)
        plt.savefig(os.path.join(output_path, con + '_heritability.png'))


def sorted_data_for_anova_analysis(input_path):
    files = sorted([i for i in os.listdir(input_path) if '.csv' in i])
    contrasts = [[i, j] for i in range(0, 3) for j in range(i + 1, i + 4) if j < 4]

    output_path = input_path + 'for_anova/'
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    for contrast in contrasts:
        data1 = pd.read_csv(input_path + files[contrast[0]], index_col=0)
        data2 = pd.read_csv(input_path + files[contrast[1]], index_col=0)
        for_anova = pd.concat((data1, data2))
        condition1 = files[contrast[0]].split('_')[0]
        condition2 = files[contrast[1]].split('_')[0]

        for_anova.to_csv(output_path + condition1 + '_vs_' + condition2 + '_pattern_similarity_for_anova.csv')


def sort_data_for_spss(input_path):
    files = sorted([i for i in os.listdir(input_path) if '.csv' in i])
    for_spss = pd.DataFrame()
    for f in files:
        data = pd.read_csv(os.path.join(input_path, f), index_col=0)
        data.columns = [i + f.split('_')[0].upper() for i in data.columns[:-2]] + ['twin', 'condition']
        for_spss = pd.concat((for_spss, data.iloc[:, :-2]), axis=1)
    for_spss['twin'] = data['twin']
    output_path = os.path.join(input_path, 'for_spss')
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    for_spss.to_csv(os.path.join(output_path, 'for_spss.csv'))


def sorted_data_for_anova_plotting(input_path):
    print('--------------')
    files = sorted([i for i in os.listdir(input_path) if '.csv' in i])
    output_path = input_path + 'for_anova_plot/'

    if not os.path.exists(output_path):
        os.mkdir(output_path)

    for_plot = pd.DataFrame()
    for f in files:
        data = pd.read_csv(input_path + f, index_col=0)
        for col in data.columns[:-2]:
            tmp = pd.DataFrame()
            tmp[['similarity', 'twin', 'condition']] = data[[col, 'twin', 'condition']]
            tmp['roi'] = [col for i in range(0, 100)]
            for_plot = pd.concat((for_plot, tmp))

    for_plot.to_csv(output_path + 'pattern_similarity_for_anova_plot.csv')
    print(output_path + 'pattern_similarity_for_anova_plot.csv')


def anova_plot(input_path, conditions, selective_dic):
    f = os.listdir(input_path)[0]
    data = pd.read_csv(input_path + f, index_col=0)
    output_path = input_path + 'anova_plot/'

    if not os.path.exists(output_path):
        os.mkdir(output_path)


    def plotting(x, y, hue, col, rois, figname):
        g = sns.factorplot(x='twin', y='similarity', hue='condition', col='roi', data=data[data.roi.isin(rois)])
        g.set_xticklabels(['mz', 'dz'], fontsize=18)
        g.set_yticklabels(fontsize=18)
        g.set_xlabels(fontsize=18)
        g.set_ylabels(fontsize=18)
        plt.savefig(output_path + figname + '.png')
        print('successfully save ' + figname + '.png')

    plotting('twin', 'similarity', 'condition', 'roi', conditions, 'system')
    for con in conditions:
        plotting('twin', 'similarity', 'condition', 'roi', selective_dic[con], con)