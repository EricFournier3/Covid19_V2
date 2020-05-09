"""
Eric Fournier 2020-04-14

"""

import os
import pandas as pd
import logging
import argparse

parser = argparse.ArgumentParser(description='Create lat_longs.tsv and ordering.tsv')
parser.add_argument('--lat-long-out', '-l', required=True, help="path to lat_longs.tsv output")
parser.add_argument('--ordering-out', '-o', required=True, help="path to ordering.tsv output")
parser.add_argument('--base-dir', '-b', required=True, help="path to working directory")
args=parser.parse_args()

logging.basicConfig(level = logging.DEBUG)

base_dir = args.base_dir

out_all_lat_long_file = args.lat_long_out
with open(out_all_lat_long_file,'w') as ll:
    pass

out_ordering = args.ordering_out
with open(out_ordering,'w') as ordering:
    pass

rta_lat_long_file = os.path.join(base_dir,"config/rta_lat_long.tsv")
rss_lat_long_file = os.path.join(base_dir,"config/rss_lat_long.tsv")
province_lat_long_file = os.path.join(base_dir,"config/province_lat_long.tsv")
canada_lat_long_file = os.path.join(base_dir,"config/canada_lat_long.tsv")
division_lat_long_file = os.path.join(base_dir,"config/division_lat_long.tsv")
country_lat_long_file = os.path.join(base_dir,"config/country_lat_long.tsv")


lat_long_level = {'rta' : ['rta','country'], 'rss' : ['rss','country'], 'division' : ['division','country'], 'country' : ['country']}
lat_long_file = {'rta' : [rta_lat_long_file], 'rss' : [rss_lat_long_file], 'division' : [division_lat_long_file,province_lat_long_file], 'country' : [country_lat_long_file,canada_lat_long_file,division_lat_long_file]} 


def Create_RTA_LatLong():
    canada_rta_file = os.path.join(base_dir,"config/canada_rta.tsv")
    canada_rta_df = pd.read_csv(canada_rta_file,delimiter="\t",header=None)
    canada_rta_df.columns = ['COUNTRY_CODE','RTA','PLACE_NAME','ADMIN_NAME_1','ADMIN_CODE_1','ADMIN_NAME_2','ADMIN_CODE_2','ADMIN_NAME_3','ADMIN_CODE_3','LATITUDE','LONGITUDE','ACCURACY']
    quebec_rta_df = canada_rta_df.loc[canada_rta_df['ADMIN_CODE_1'] == 'QC',['RTA','LATITUDE','LONGITUDE']]

    quebec_rta_df.to_csv(rta_lat_long_file,sep="\t",index=False,header=None)
         

Create_RTA_LatLong()

for level, sublevel in lat_long_level.items():
    print("************************* In Level " + level)
    for sub in sublevel:
        for myfile in lat_long_file[sub]:
            df_lat_long = None
            #print("MY FILE IS ",myfile)
            df_lat_long = pd.read_csv(myfile,delimiter="\t",header=None)
            df_lat_long.insert(loc=0,column="",value=level,allow_duplicates=True)
            #print(df_lat_long)
            df_ordering = df_lat_long.iloc[:,0:2]
            df_lat_long.to_csv(out_all_lat_long_file,sep="\t",index=False,header=False,mode='a')
            df_ordering.to_csv(out_ordering,sep="\t",index=False,header=False,mode='a')
            

logging.info("End")
