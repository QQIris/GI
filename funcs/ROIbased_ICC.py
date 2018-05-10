import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import data_preprocessing


def gss_get_avg_data(input_path, gss_mask_input, conditions, vtx_map, twinid, output_path):
    print('get average data under gss mask')
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    beta_files = [i for i in os.listdir(input_path) if '.csv' in i]
    gss_files = sorted([i for i in os.listdir(gss_mask_input) if '.csv' in i])
    gss_file_set = [[gss_files[i],gss_files[i+len(gss_files)/2]] for i in range(0,len(gss_files)/2)]
    rois = [i.split('_')[1] for i in gss_files[:len(gss_files)/2]]
    for con in conditions:
        f = [i for i in beta_files if con in i][0]
        print('-------------------------------------')
        print(f)
        data = pd.read_csv(input_path+f, index_col=0)
        data_l = data.iloc[:,:29696]
        data_r = data.iloc[:,29696:59412]
        data_l.columns = vtx_map[0]
        data_r.columns = vtx_map[1]
        avg_data = pd.DataFrame(index=data.index, columns=rois)
        for sets in gss_file_set:
            gss_mask_l = pd.read_csv(os.path.join(gss_mask_input,sets[0]), index_col=0)
            gss_mask_r = pd.read_csv(os.path.join(gss_mask_input,sets[1]), index_col=0)
            roi = sets[0].split('_')[1]
            print(roi)
            for family in twinid.index:
                print(family)
                avg_l = data_l.ix[list(twinid.ix[family])][gss_mask_l.ix[family].dropna()]
                avg_r = data_r.ix[list(twinid.ix[family])][gss_mask_r.ix[family].dropna()]
                avg = pd.concat((avg_l, avg_r),axis=1).mean(axis=1)
                avg_data[roi].loc[avg.index] = avg

        avg_data.to_csv(os.path.join(output_path , '_'.join(f.split('_')[:-1]) + '_gss_roi_mean.csv'))
        print('successful save csv file')
    return  output_path

def get_avg_data(input_path, vtx_mask, vtx_map, all_regions, output_path):
    files = sorted([i for i in os.listdir(input_path) if '_fix_' in i])
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    for f in files:
        print(f)
        con = f.split('_')[0]
        data = pd.read_csv(os.path.join(input_path, f), index_col=0)
        data_l = data.iloc[:,:29696]
        data_r = data.iloc[:,29696:59412]
        data_l.columns = vtx_map[0]
        data_r.columns = vtx_map[1]
        avg_data = pd.DataFrame(index=data.index)
        for region in all_regions:
            masked_data_l = data_l[vtx_mask[0].loc[region].dropna()]
            masked_data_r = data_r[vtx_mask[1].loc[region].dropna()]
            masked_data = pd.concat((masked_data_l, masked_data_r), axis=1)
            avg_data[region] = masked_data.mean(axis=1)
        avg_data.to_csv(os.path.join(output_path, '_'.join(f.split('_')[:-1]) + '_roi_mean.csv'))
        print('successfully save data')

def get_topn_avg_data(input_path, vtx_mask, vtx_map, all_regions, output_path, npercent=0.2):
    files = sorted([i for i in os.listdir(input_path) if '_avg_' in i])
    print(npercent)
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    for f in files:
        print(f)
        data = pd.read_csv(os.path.join(input_path, f), index_col=0)
        data_l = data.iloc[:,:29696]
        data_r = data.iloc[:,29696:59412]
        data_l.columns = vtx_map[0]
        data_r.columns = vtx_map[1]
        avg_data = pd.DataFrame(index=data.index)

        for region in all_regions:
            print(region)
            mask_l = vtx_mask[0].loc[region].dropna()
            mask_r = vtx_mask[1].loc[region].dropna()

            masked_data_l = data_l[mask_l]
            masked_data_r = data_r[mask_r]
            roi_avg = []

            for sub in masked_data_l.index:
                sub_l = masked_data_l.loc[sub].nlargest(int(mask_l.count() * npercent))
                sub_r = masked_data_r.loc[sub].nlargest(int(mask_r.count() * npercent))
                roi_avg.append(np.mean(np.concatenate((sub_r, sub_r))))

            avg_data[region] = roi_avg
        avg_data.to_csv(os.path.join(output_path, '_'.join(f.split('_')[:-1]) + '_roi_topn_mean.csv'))
        print('successfully save data')


def rm_outliers_avg(input_path, mzid, dzid, thr=[-2,2]):
    files = sorted([i for i in os.listdir(input_path) if '.csv' in i])
    output_path = os.path.join(input_path, 'rm_outliers', 'thr'+str(thr[1]))
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    for f in files:
        data = pd.read_csv(os.path.join(input_path, f), index_col=0)
        mzid_list = np.concatenate((mzid['twin1'].values, mzid['twin2'].values))
        dzid_list = np.concatenate((dzid['twin1'].values, dzid['twin2'].values))
        mzdata_ready = data_preprocessing.data_preparation(data.loc[mzid_list])
        dzdata_ready = data_preprocessing.data_preparation(data.loc[dzid_list])
        mzdata_ready.all_roi_remove_outleirs(thr)
        dzdata_ready.all_roi_remove_outleirs(thr)
        mzdata_ready.rm_twin_pair(mzid)
        dzdata_ready.rm_twin_pair(dzid)

        rm_data = pd.concat((mzdata_ready.data, dzdata_ready.data))

        rm_data.to_csv(os.path.join(output_path, f[:-4]+'_rm_outliers.csv'))

def gss_avg_data_distribution(input_path, selective_dic, colors_dic):
    print('------------------------------------------')
    print('generate hist of averaged data')

    files = sorted([i for i in os.listdir(input_path) if '.csv' in i])
    output_path = os.path.join(input_path, 'hist/')
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    for f in files:
        con = f.split('_')[0]
        regions = [con] + selective_dic[con]
        data = pd.read_csv(os.path.join(input_path, f), index_col=0)
        data[regions].plot(kind='hist', subplots=True, layout=(1, len(regions)), figsize=(len(regions) * 5, 5),
                           sharey=True, color=colors_dic[con])
        plt.subplots_adjust(hspace=0.10, wspace=0.10)
        plt.savefig(output_path + f[:-4] + '_hist.png')
        plt.close()

def sorted_avg_data_in_twin_forms(input_path, mzid, dzid):
    print('------------------------------------------------------')
    print('we are now sorted averaged data into twin forms')
    files = sorted([i for i in os.listdir(input_path) if '.csv' in i])
    output_path = os.path.join(input_path, 'twin_form/')
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    for f in files:
        print(f)
        data = pd.read_csv(os.path.join(input_path + f), index_col=0)
        data_sort_mz = data_preprocessing.twin_data_preparation(mzid, in_twin_forms=True)
        data_sort_dz = data_preprocessing.twin_data_preparation(dzid, in_twin_forms=True)
        for col in data.columns:
            data_sort_mz.load_univariate_twin_data(data[col], col, merge_to_twindata=True)
            data_sort_dz.load_univariate_twin_data(data[col], col, merge_to_twindata=True)

        data_sort_mz.data.to_csv(os.path.join(output_path, f.split('_')[0] + '_condition_sorted_data_into_mztwin.csv'))
        data_sort_dz.data.to_csv(os.path.join(output_path, f.split('_')[0] + '_condition_sorted_data_into_dztwin.csv'))


    return output_path

def sig_icc_results(input_path):
    files = sorted([i for i in os.listdir(input_path) if '.csv' in i])
    for n, f in enumerate(files):
        data = pd.read_csv(input_path + f, index_col=0)
        con = f.split('_')[0]
        data['mz-dz'] = data['iccmz'] - data['iccdz']
        data.to_csv(os.path.join(input_path, f))

def des_icc(input_path, conditions, selective_dic):

    print('-----------------------------------')
    print('we are now display describtive results of icc')

    files = sorted([i for i in os.listdir(input_path) if '.csv' in i])
    fig, axes = plt.subplots(4, 5, figsize=(20, 20), sharex='col', sharey='all')
    colors = [['b','b','b','b'],
              ['g','g'],
              ['r','r','r'],
              ['c','c','c']]
    data_con = pd.DataFrame()

    for n, f in enumerate(files):

        data = pd.read_csv(input_path+f, index_col=0)
        con = f.split('_')[0]

        data['iccmz'].ix[selective_dic[con]].plot(kind='bar', ax=axes[0,n+1], color=colors[n], fontsize=18,
                                                  grid=True, width=1)
        data['iccdz'].ix[selective_dic[con]].plot(kind='bar', ax=axes[1,n+1], color=colors[n], fontsize=18,
                                                  grid=True, width=1)
        data['mz-dz'].ix[selective_dic[con]].plot(kind='bar', ax=axes[2,n+1], color=colors[n], fontsize=18,
                                                  grid=True, width=1)
        data['sigmz-dz'].ix[selective_dic[con]].plot(kind='bar', ax=axes[3,n+1], color=colors[n], fontsize=18,
                                                  grid=True, width=1)

        data_con[con] = data[['iccmz','iccdz','mz-dz','sigmz-dz']].ix[con]

    data_con = data_con.T
    data_con['iccmz'].ix[conditions].plot(kind='bar', ax=axes[0,0], color=['b','g','r','c'], grid=True, width=1, fontsize=18,)
    data_con['iccdz'].ix[conditions].plot(kind='bar', ax=axes[1,0], color=['b','g','r','c'], grid=True, width=1, fontsize=18,)
    data_con['mz-dz'].ix[conditions].plot(kind='bar', ax=axes[2,0], color=['b','g','r','c'], grid=True, width=1, fontsize=18,)
    data_con['sigmz-dz'].ix[conditions].plot(kind='bar', ax=axes[3,0], color=['b','g','r','c'], grid=True, width=1, fontsize=18,)

    fig.subplots_adjust(hspace=0.11, wspace=0, top=0.95, bottom=0.15)
    output_path = os.path.join(input_path, 'bar/')
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    # [ax.set_xlabel(fontsize=18) for ax in axes.ravel()]
    # [ax.set_ylabel(fontsize=18) for ax in axes.ravel()]
    plt.subplots_adjust(wspace=0)
    plt.savefig(os.path.join(output_path,'icc_results_desribe.png'))
    plt.close()

def des_cate_icc(input_path, conditions, selective_dic):

    print('-----------------------------------')
    print('we are now display category results of icc')
    files = sorted([i for i in os.listdir(input_path) if '.csv' in i])
    colors = ['b','g','r','c']
    data_con = pd.DataFrame()

    output_path = os.path.join(input_path, 'bar_cate/')
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    fig1, axes1 = plt.subplots(4, 4, figsize=(20, 20), sharex='col', sharey='all')
    for r,region in enumerate(conditions):

        all_sys_data = pd.DataFrame()
        for n, con in enumerate(conditions):
            f = [i for i in files if con in i][0]
            data = pd.read_csv(input_path + f, index_col=0)
            all_sys_data[con] = data[['iccmz', 'iccdz', 'mz-dz', 'sigmz-dz', 'sigmz-posdz']].ix[region]

        all_sys_data = all_sys_data.T
        all_sys_data['iccmz'].plot(kind='bar', ax=axes1[0, r], color=colors, grid=True, width=1, fontsize=18)
        all_sys_data['iccdz'].plot(kind='bar', ax=axes1[1, r], color=colors, grid=True, width=1, fontsize=18)
        all_sys_data['mz-dz'].plot(kind='bar', ax=axes1[2, r], color=colors, grid=True, width=1, fontsize=18)
        all_sys_data['sigmz-dz'].plot(kind='bar', ax=axes1[3, r], color=colors, grid=True, width=1, fontsize=18)
    fig1.subplots_adjust(wspace=0)
    fig1.savefig(os.path.join(output_path, 'all_system_category_results.png'))

class icc_display():

    def __init__(self, input_path, conditons, colors, selective_dic):
        self.input_path = input_path
        self.files = sorted([i for i in os.listdir(self.input_path) if '.csv' in i])
        self.conditons = conditons
        self.selective_dic = selective_dic
        self.output_path = os.path.join(input_path, 'bar/')
        self.colors = colors

    def describe_icc_twin_seperate(self):
        print('-----------------------------------')
        print('we are now display describtive results of icc')

        files = sorted([i for i in os.listdir(self.input_path) if '.csv' in i])
        fig, axes = plt.subplots(4, 5, figsize=(20, 20), sharex='col', sharey='all')

        data_con = pd.DataFrame()

        for n, f in enumerate(files):
            data = pd.read_csv(self.input_path + f, index_col=0)
            con = f.split('_')[0]

            data['iccmz'].ix[self.selective_dic[con]].plot(kind='bar', ax=axes[0, n + 1], color=self.colors[n], fontsize=18,
                                                      grid=True, width=1)
            data['iccdz'].ix[self.selective_dic[con]].plot(kind='bar', ax=axes[1, n + 1], color=self.colors[n], fontsize=18,
                                                      grid=True, width=1)
            data['mz-dz'].ix[self.selective_dic[con]].plot(kind='bar', ax=axes[2, n + 1], color=self.colors[n], fontsize=18,
                                                      grid=True, width=1)
            data['sigmz-dz'].ix[self.selective_dic[con]].plot(kind='bar', ax=axes[3, n + 1], color=self.colors[n], fontsize=18,
                                                         grid=True, width=1)

            data_con[con] = data[['iccmz', 'iccdz', 'mz-dz', 'sigmz-dz']].ix[con]

        data_con = data_con.T
        data_con['iccmz'].ix[self.conditions].plot(kind='bar', ax=axes[0, 0], color=['b', 'g', 'r', 'c'], grid=True, width=1,
                                              fontsize=18, )
        data_con['iccdz'].ix[self.conditions].plot(kind='bar', ax=axes[1, 0], color=['b', 'g', 'r', 'c'], grid=True, width=1,
                                              fontsize=18, )
        data_con['mz-dz'].ix[self.conditions].plot(kind='bar', ax=axes[2, 0], color=['b', 'g', 'r', 'c'], grid=True, width=1,
                                              fontsize=18, )
        data_con['sigmz-dz'].ix[self.conditions].plot(kind='bar', ax=axes[3, 0], color=['b', 'g', 'r', 'c'], grid=True,
                                                 width=1, fontsize=18, )

        fig.subplots_adjust(hspace=0.11, wspace=0, top=0.95, bottom=0.15)
        output_path = os.path.join(self.input_path, 'bar/')
        if not os.path.exists(output_path):
            os.mkdir(output_path)

        plt.subplots_adjust(wspace=0)
        plt.savefig(os.path.join(output_path, 'icc_results_desribe.png'))
        plt.close()

    def describe_icc_twin(self,cons):
        icc_results = pd.DataFrame()

        for con in cons:
            f = [i for i in self.files if con in self.files][0]
            data = pd.read_csv(os.path.join(self.input_path, f), index_col=0)
            icc_results[con] = data[['iccmz','iccdz']].loc[con]
        icc_results.T.plot.bar(grid=True, width=0.65, color=['black', 'gray'])
        if os.path.exists(self.output_path):
            os.mkdir(self.output_path)
        plt.savefig(os.path.join(self.output_path,'icc_results_twin_system.png'))


    def cate_icc(self, regions):
        pass
