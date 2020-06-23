
import os
import sys
import csv
import pandas as pd
import numpy as np
import utils_full_network as utils

data_loc = "6D_UC"

output_file_name = data_loc + "_output.csv"
output_file_name


columns = ['cap_reserve_margin','voll','voso','NG_price','Wind price','line_ex','FS_cost','SS_cost','max_wind','avg_wind','max_gas','filename']
output_agg_df = pd.DataFrame(index=range(1,len(os.listdir(data_loc))+1),columns=columns)
# load stacks_df
stacks_df = pd.read_csv('./system_data/stacks_df_full_year.csv',header=0,index_col=0,parse_dates=True)
# load system_gens
system_gens_raw = pd.read_csv('./system_data/system_gens.csv',header=0,index_col=0,parse_dates=True)
system_gens_orginal = pd.read_csv('./RTS_Data/SourceData/gen.csv',header=0,index_col=0,parse_dates=True)
system_gens, stacks_df  = utils.green_field_rts(system_gens_raw,stacks_df, extra_wind=True)
system_gens = utils.reset_all_gens(system_gens)
system_gens = utils.add_RTS_thermal_only(system_gens, system_gens_orginal)
system_gens= utils.add_wind_sites(system_gens)
system_gens = utils.retire_units(system_gens)
system_gens = utils.add_new_thermal(system_gens, system_gens_orginal)

def post_process_builds(system_gens,filename):
    x = []
    solution_df = pd.read_csv(filename, header=None)
    solution_df.columns = ['stage', 'node', 'name', 'index', 'value']
    solution_df.value.fillna(0, inplace=True)
    df_fs=solution_df.query("stage=='FirstStage'").copy()
    a_list = []
    for item in df_fs['index'].values:
        a_list.append(item.strip())

    system_gens = system_gens.loc[a_list,:]


    df_build=pd.DataFrame(columns=['number on system','new number on system'],index=system_gens.index)
    df_build['number on system']=system_gens['number on system']
    list_a=[]
    for i in df_fs['value'].astype(float).values:
        list_a.append(i)
    df_build['new number on system']= list_a
    df_build['new builds']=list_a
    df_build['category']=system_gens['gen cat']
    df_build['max gen']=system_gens['max gen']
    df_build['average available capacity']=system_gens['max gen']*system_gens['capacity discount']*df_build['new builds']
    df_build['max available capacity']=system_gens['max gen']*df_build['new builds']
    df_build['max build available capacity']=system_gens['max gen']*system_gens['capacity discount']
    df_build['busID']=system_gens['busID']

    df_build = df_build.loc[df_build['new number on system'] == 1.0]

    wind_builds = df_build.loc[df_build['category']=="Wind"]

    wind_max_gen = sum(list(wind_builds['max gen']))
    x.append(wind_max_gen)
    wind_avg_gen = sum(list(wind_builds['average available capacity']))
    x.append(wind_avg_gen)

    gas_builds =  df_build.loc[df_build['category'].str.contains('Gas')]

    gas_max_gen = sum(list(gas_builds['max gen']))
    x.append(gas_max_gen)
    return x

# output_agg_df.iloc[0,0:4] = x
# output_agg_df.iloc[0,4] = 1
# output_agg_df

# # Data frames for aggregation
# output_df = pd.read_csv("test_run/results_35/"+"ph_StageCostDetail.csv")
# input_df = pd.read_csv("test_run/results_35/"+"parameters_df_sample.csv")
#
# x = []
# # process inputs
# num_scenarions = input_df.loc[0, 'number of scenarios']
# cap_reserve_margin = input_df.loc[0, 'capacity reserves']
# voll = input_df.loc[0, 'loss of load cost']
# voso = input_df.loc[0, 'system overload cost']
# NG_price = input_df.loc[0, 'NG price']
# # scale back to [-1,1]
# cap_reserve_margin = 2*(cap_reserve_margin-0.1)/(0.2-0.1)-1
# voll = 2*(voll - 5000)/(15000-5000)-1
# voso = 2*(voso-100)/(1000-100)-1
# NG_price = 2*(NG_price-2.5)/(4.7-2.5)-1
# x.append(cap_reserve_margin)
# x.append(voll)
# x.append(voso)
# x.append(NG_price)
# output_agg_df.iloc[0,0:4] = x
#
# #Process outputs
# FS_cost = output_df.iloc[0,-1]
# output_agg_df.iloc[file_index,4] = FS_cost
# my_column = list(output_df.iloc[num_scenarions-1:2*num_scenarions-1, -1])
# SS_cost = sum(my_column)
# output_agg_df.iloc[file_index,5] = SS_cost
file_index = 0
for filename in os.listdir(data_loc):
    if not filename.startswith('.'):
        print filename
        output_df = pd.read_csv(data_loc+'/'+filename+'/'+"ph_StageCostDetail.csv")
        input_df = pd.read_csv(data_loc+'/'+filename+'/'+"parameters_df_sample.csv")
        x = []
        # process inputs
        num_scenarions = input_df.loc[0, 'number of scenarios']
        cap_reserve_margin = input_df.loc[0, 'capacity reserves']
        voll = input_df.loc[0, 'loss of load cost']
        voso = input_df.loc[0, 'system overload cost']
        NG_price = input_df.loc[0, 'NG price']
        Wind_price = input_df.loc[0, 'Wind price']
        line_ex = input_df.loc[0, 'line_ex']
        weight = input_df.loc[0, 'weight']
        # scale back to [-1,1]
        cap_reserve_margin = 2*(cap_reserve_margin-0.1)/(0.2-0.1)-1
        voll = 2*(voll - 4000)/(15000-4000)-1
        voso = 2*(voso-100)/(1000-100)-1
        NG_price = 2*(NG_price-1.6)/(8-1.6)-1
        Wind_price = 2*(Wind_price - 0.85)/(1.15-0.85)-1
        line_ex = 2*(line_ex-0.5)/(1.5-0.5)-1
        x.append(cap_reserve_margin)
        x.append(voll)
        x.append(voso)
        x.append(NG_price)
        x.append(Wind_price)
        x.append(line_ex)
        # Save inputs to data frame
        output_agg_df.iloc[file_index, 0:6] = x


        # Process outputs
        FS_cost = output_df.iloc[0, -1]
        my_column = list(output_df.iloc[num_scenarions-1:2*num_scenarions-1, -1])
        SS_cost = sum(my_column)*365/num_scenarions/weight
        post_process_filename = data_loc+'/'+filename+'/'+"ph.csv"
        xx = post_process_builds(system_gens,post_process_filename)

        # Save outputs to data frame
        output_agg_df.iloc[file_index, 6] = FS_cost
        output_agg_df.iloc[file_index, 7] = SS_cost
        output_agg_df.iloc[file_index, 8:11] = xx
        output_agg_df.iloc[file_index, 11] = filename
        file_index = file_index + 1
        print file_index

output_agg_df
output_agg_df.to_csv(output_file_name)
