import pandas as pd
import os

#pour canada_rta.tsv voir http://download.geonames.org/export/zip/

pd.set_option('display.max_rows', None)

base_dir = "/data/Applications/GitScript/Covid19_V2/NextStrainFiles/config"

lspq_sgil_extract = "/data/PROJETS/COVID_19/TestOnlyQuebec2/lspq/sgil_extract.tsv" 
df_lspq = pd.read_csv(lspq_sgil_extract,delimiter="\t",index_col=False)
rta_in_df_lspq = df_lspq['POSTAL_CODE'].unique()
rta_in_df_lspq = [x.lower() for x in rta_in_df_lspq]
#print(rta_in_df_lspq)

canada_rta_file = os.path.join(base_dir,"canada_rta.tsv")

canada_rta_df = pd.read_csv(canada_rta_file,delimiter="\t",header=None)
canada_rta_df.columns = ['COUNTRY_CODE','RTA','PLACE_NAME','ADMIN_NAME_1','ADMIN_CODE_1','ADMIN_NAME_2','ADMIN_CODE_2','ADMIN_NAME_3','ADMIN_CODE_3','LATITUDE','LONGITUDE','ACCURACY']
quebec_rta_df = canada_rta_df.loc[canada_rta_df['ADMIN_CODE_1'] == 'QC',['RTA','LATITUDE','LONGITUDE']]
#print(len(quebec_rta_df.index))
quebec_rta_values = quebec_rta_df.iloc[:,0].values
quebec_rta_values = [x.lower() for x in quebec_rta_values]
#print(quebec_rta_values)

missing_rta = [rta for rta in rta_in_df_lspq if rta not in quebec_rta_values]

#print(missing_rta)
#['a9b', 'j3j', 'l6n', 'j3s', 'j9w', 'k4k', 'j6l', 'j5s', 'k1k', 'k1h', 'j8j', 'k1c', 'k1l', 'k2c']
#print(df_lspq[(~df_lspq['POSTAL_CODE'].isin([x.upper() for x in missing_rta])) & (~df_lspq['NO_LSPQ'].isin(['L00256307']))])
#print(df_lspq[df_lspq['POSTAL_CODE'].isin(['G1H','J2H'])])
print(df_lspq[df_lspq['POSTAL_CODE'].isin([])])

#L00256307
