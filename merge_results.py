import pandas as pd
import numpy as np 
import matplotlib as mpl
import matplotlib.pyplot as plt
import pickle
import glob

### Loading in all Data
# interface_df = pd.read_csv("/qfs/projects/proteometer/pheno_analysis/human_proteome_precaculated/interfaces_precaculated.csv")
# pocket_df = pd.read_csv("/qfs/projects/proteometer/pheno_analysis/human_proteome_precaculated/pockets_precaculated.csv")
sasa_data = "/qfs/projects/proteometer/pheno_analysis/sasa_csv"
charge_data = "/qfs/projects/proteometer/pheno_analysis/charge_csv"
# stability_data = pd.read_csv("/rcfs/projects/proteometer/human_proteome_precaculated/stability_precaculated.csv")
# pka_data = pd.read_csv("/rcfs/projects/proteometer/human_proteome_precaculated/pka_precaculated.csv")


sasa_files = glob.glob(sasa_data+ "/*")
charge_files = glob.glob(charge_data+ "/*")


# merge all pocket files
sasa_list = []
for file in sasa_files:
    sasa_list.append(pd.read_csv(file))
sasa_df = pd.concat(sasa_list)

sasa_df.to_csv("/rcfs/projects/proteometer/human_proteome_precaculated/sasa_precaculated.csv")

# merge all charge files
charge_list = []
for file in charge_files:
    charge_list.append(pd.read_csv(file))
charge_df = pd.concat(charge_list)

charge_df.to_csv("/rcfs/projects/proteometer/human_proteome_precaculated/charge_precaculated.csv")


# sasa_charge =  pd.merge(sasa_df, charge_df[['charge', 'charge_method', 'UNIPROT', 'RES_NUM']], how= 'inner', left_on=['UNIPROT', 'RES_NUM'], right_on=['UNIPROT', 'RES_NUM'])
# sasa_charge_pockets = pd.merge(sasa_charge, pocket_df[['protein_acc', 'wild_type', 'position', 'plddt','foldx_ddg_min','foldx_ddg_max','foldx_ddg_abs_max','foldx_ddg_abs_median','closest_pocket','inside_pocket', 'min_distance_from_pocket','mean_distance_from_pocket']], how= 'inner', left_on=['UNIPROT', 'RES_NUM'], right_on=['protein_acc', 'position'])
# sasa_charge_pockets_interfaces =  pd.merge(sasa_charge_pockets, interface_df[['protein_acc', 'position', 'closest_interface','inside_interface', 'min_distance_from_interface','mean_distance_from_interface']], how= 'inner', left_on=['UNIPROT', 'RES_NUM'], right_on=['protein_acc', 'position'])


# sasa_charge_pockets_interfaces.to_csv("/rcfs/projects/proteometer/human_proteome_precaculated/all_combined/all_combined.csv")

