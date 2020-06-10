import os
from datetime import date
import pandas as pd

today = date.today()
today = today.strftime("%Y%m%d")

out_res_dir = "/data/PROJETS/COVID_19/Missing_SGIL_RTA/"
in_nextstrain_dir = "/data/Applications/GitScript/Covid19_V2/NextStrainFiles/"

missing_sgil_rta_out_file="missing_sgil_rta_" + today + ".xls"
missing_sgil_rta_out_path=os.path.join(out_res_dir,missing_sgil_rta_out_file)
lspq_sgil_extract = os.path.join(in_nextstrain_dir,"data/lspq","sgil_extract.tsv")
canada_rta_file = os.path.join(in_nextstrain_dir,"config","canada_rta.tsv")

rta_lat_long_file = os.path.join(out_res_dir,"rta_lat_long.tsv")

out_lspq_partage = os.path.join("/mnt/Partage/LSPQ_Partage/DiagnosticMoleculaire/PROJETS/Covid19/MissingRTA")

def Create_RTA_LatLong():
    canada_rta_df = pd.read_csv(canada_rta_file,delimiter="\t",header=None)
    canada_rta_df.columns = ['COUNTRY_CODE','RTA','PLACE_NAME','ADMIN_NAME_1','ADMIN_CODE_1','ADMIN_NAME_2','ADMIN_CODE_2','ADMIN_NAME_3','ADMIN_CODE_3','LATITUDE','LONGITUDE','ACCURACY']
    quebec_rta_df = canada_rta_df.loc[canada_rta_df['ADMIN_NAME_1'].isin(['Ontario','Quebec'])  ,['RTA','LATITUDE','LONGITUDE']]

    quebec_rta_df.to_csv(rta_lat_long_file,sep="\t",index=False,header=None)


def CheckMissingRTA():
    df_lspq = pd.read_csv(lspq_sgil_extract,delimiter="\t",index_col=False)
    rta_df = pd.read_csv(rta_lat_long_file,delimiter="\t",header=None)
    rta_val = rta_df.iloc[:,0].values
    rta_val = [x.lower() for x  in rta_val]

    rta_in_df_lspq = df_lspq['POSTAL_CODE'].unique()
    rta_in_df_lspq = [x.lower() for x in rta_in_df_lspq]

    missing_rta_in_latlong = [rta for rta in rta_in_df_lspq if rta not in rta_val]
    missing_rta_in_latlong = [x.upper() for x in missing_rta_in_latlong]

    missing_rta = df_lspq.loc[df_lspq['POSTAL_CODE'].isin(missing_rta_in_latlong),['NO_LSPQ','POSTAL_CODE','RSS_PATIENT']]
    print(missing_rta)
    missing_rta.to_excel(missing_sgil_rta_out_path,index=False,sheet_name='Missing rta')
    os.system("sudo cp " + missing_sgil_rta_out_path + " " + out_lspq_partage)



Create_RTA_LatLong()
CheckMissingRTA()


print("End")
