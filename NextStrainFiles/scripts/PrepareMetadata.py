"""
Eric Fournier 2020-04-26

"""

import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
import re
import logging
from googletrans import Translator
import urllib3
import argparse


parser = argparse.ArgumentParser(description='Create metadata.tsv and sequences.fasta')
parser.add_argument('--base-dir', '-b', required=True, help="path to working directory")
parser.add_argument('--metadata-out', '-m', required=True, help="path to resulting metadata.tsv")
parser.add_argument('--sequences-out', '-s', required=True, help="path to resulting sequences.fasta")
args=parser.parse_args()


logging.getLogger("urllib3").setLevel(logging.WARNING)
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

pd.set_option('display.max_columns', 20)

translator = Translator()
logging.basicConfig(level = logging.DEBUG)

base_dir = args.base_dir
lat_long_file = os.path.join(base_dir,"config/lat_longs.tsv")
country_lat_long_file = os.path.join(base_dir,"config/country_lat_long.tsv")
ordering_file = os.path.join(base_dir,"config/ordering.tsv")

gisaid_metadata = os.path.join(base_dir,"gisaid/metadata.tsv") 
lspq_sgil_extract = os.path.join(base_dir,"lspq/sgil_extract.tsv")

gisaid_ref_sequences = os.path.join(base_dir,"gisaid/gisaid_wuhan_ref_20200425.fasta")
gisaid_sequences = os.path.join(base_dir,"gisaid/gisaid_all.fasta")
lspq_sequences = os.path.join(base_dir,"lspq/sequences.fasta")

df_lspq = pd.read_csv(lspq_sgil_extract,delimiter="\t",index_col=False)
df_gisaid = pd.read_csv(gisaid_metadata,delimiter="\t",index_col=False)
df_lat_long = pd.read_csv(lat_long_file,delimiter='\t',header=None)

rec_list = []
rec_id_list = []

out_metadata = args.metadata_out
out_sequences = args.sequences_out


def ConvertFrench2English(countries_list):

    still_missing_countries = []


    for country in countries_list:
        english_country = str(translator.translate(country,dest='en').text).capitalize()
        logging.info("Convert " + country + " in english => " + english_country)
        if len(df_lat_long.loc[df_lat_long[1] == english_country,1].values) > 0:
            df_lspq.loc[df_lspq['VOYAGE_PAYS_1'].isin([str(country).capitalize(),str(country).lower(),str(country).upper()]),'VOYAGE_PAYS_1'] = english_country
        else:
            still_missing_countries.append(country)
    return(still_missing_countries.copy())


def CheckMissingCountry():
    locs_w_latlong_orig = df_lat_long.loc[df_lat_long[0]=='country',1].values 
    locs_w_latlong = [x.lower() for x in locs_w_latlong_orig]

    locs_in_lspq_orig = [x for x in df_lspq['VOYAGE_PAYS_1'].unique() if x != "AUCUN_VOYAGE"]
    locs_in_lspq = [x.lower() for x in locs_in_lspq_orig]

    missing_latlong_locs = [loc for loc in locs_in_lspq if loc not in locs_w_latlong]
    still_missing_country_in_latlong = ConvertFrench2English(missing_latlong_locs)

    if len(still_missing_country_in_latlong) > 0:
        return still_missing_country_in_latlong
    else:
        return None

missing_country = CheckMissingCountry()

if  missing_country:
    logging.error(" ----------- Countries missing -------------")
    logging.info("Add the following countries (" + str(len(missing_country)) + ") in "  + "/data/Applications/GitScript/Covid19/NextStrainFiles/config/MissingSgilCountries.tsv")
    logging.info(str(missing_country))
    logging.info(" and run again")
    exit(1)    

id_ref_list = []

for rec in SeqIO.parse(gisaid_ref_sequences,'fasta'):
    try:
        seqid = re.search(r'^\S+\/(\S+\/\S+\/\d{4})\|\S+',rec.id).group(1)
    except:
        seqid = re.search(r'^\S+\/(\S+\/\d{4})\|\S+',rec.id).group(1)

    id_ref_list.append(seqid)
    rec_id_list.append(seqid)
    rec.id = seqid
    rec.name = seqid
    rec.description = ""
    rec_list.append(rec)

for rec in SeqIO.parse(gisaid_sequences,'fasta'):
    
    if rec.id not in id_ref_list:
        try:
            seqid = re.search(r'^(\S+\/\S+\/\d{4})',rec.id).group(1)
        except:
            print("BUG WITH " + rec.id)
            seqid = re.search(r'^\S+\/(\S+\/\S+\/\d{4})\|\S+',rec.id).group(1)

        rec_id_list.append(seqid)
        rec.id = seqid
        rec.name = seqid
        rec.description = ""
        rec_list.append(rec)

for rec in SeqIO.parse(lspq_sequences,'fasta'):
    seqid = re.search(r'(^\S+)\/\S+\/\S+',rec.id).group(1)
    rec_id_list.append(seqid)
    rec.id = seqid
    rec.name = seqid
    rec.description = ""
    rec_list.append(rec)

subset_lspq = df_lspq[df_lspq['NO_LSPQ'].isin(rec_id_list)]
subset_gisaid = df_gisaid[df_gisaid['strain'].isin(rec_id_list)]

#subset_lspq_subcol = subset_lspq[['NO_LSPQ','DATE_PRELEV','SEX','AGE','RSS_PATIENT','POSTAL_CODE','VOYAGE_PAYS_1']] this create a view https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy
subset_lspq_subcol = subset_lspq.loc[:,('NO_LSPQ','DATE_PRELEV','SEX','AGE','RSS_PATIENT','POSTAL_CODE','VOYAGE_PAYS_1','MAX_CT')] # this create a copy
subset_lspq_subcol.columns = ['strain','date','sex','age','rss','rta','voyage','ct']
subset_lspq_subcol.insert(loc=1,column='virus',value='ncov',allow_duplicates=True)
subset_lspq_subcol.insert(loc=5,column='country',value='Canada',allow_duplicates=True)
subset_lspq_subcol.insert(loc=6,column='division',value='Quebec',allow_duplicates=True)
subset_lspq_subcol.insert(loc=9,column='country_exposure',value='',allow_duplicates=True)
subset_lspq_subcol.insert(loc=10,column='division_exposure',value='',allow_duplicates=True)
subset_lspq_subcol.insert(loc=11,column='rss_exposure',value='',allow_duplicates=True)
subset_lspq_subcol.insert(loc=12,column='rta_exposure',value='',allow_duplicates=True)
subset_lspq_subcol.loc[subset_lspq_subcol.voyage == 'AUCUN_VOYAGE','country_exposure'] = subset_lspq_subcol.country
subset_lspq_subcol.loc[subset_lspq_subcol.voyage == 'AUCUN_VOYAGE','division_exposure'] = subset_lspq_subcol.division
subset_lspq_subcol.loc[subset_lspq_subcol.voyage == 'AUCUN_VOYAGE','rss_exposure'] = subset_lspq_subcol.rss
subset_lspq_subcol.loc[subset_lspq_subcol.voyage == 'AUCUN_VOYAGE','rta_exposure'] = subset_lspq_subcol.rta
subset_lspq_subcol.loc[subset_lspq_subcol.voyage != 'AUCUN_VOYAGE','country_exposure'] = subset_lspq_subcol.voyage
subset_lspq_subcol.loc[subset_lspq_subcol.voyage != 'AUCUN_VOYAGE','division_exposure'] = subset_lspq_subcol.voyage
subset_lspq_subcol.loc[subset_lspq_subcol.voyage != 'AUCUN_VOYAGE','rss_exposure'] = subset_lspq_subcol.voyage
subset_lspq_subcol.loc[subset_lspq_subcol.voyage != 'AUCUN_VOYAGE','rta_exposure'] = subset_lspq_subcol.voyage

del subset_lspq_subcol['voyage']
#print(subset_lspq_subcol)

subset_gisaid_subcol = subset_gisaid.loc[:,('strain','virus','date','sex','age','country','division','country_exposure','division_exposure')]
subset_gisaid_subcol.insert(loc=7,column='rss',value='',allow_duplicates=True)
subset_gisaid_subcol.insert(loc=8,column='rta',value='',allow_duplicates=True)
subset_gisaid_subcol.insert(loc=11,column='rss_exposure',value='',allow_duplicates=True)
subset_gisaid_subcol.insert(loc=12,column='rta_exposure',value='',allow_duplicates=True)
subset_gisaid_subcol.insert(loc=13,column='ct',value='0',allow_duplicates=True)
subset_gisaid_subcol.loc[:,'rss'] = subset_gisaid_subcol.country
subset_gisaid_subcol.loc[:,'rta'] = subset_gisaid_subcol.country
subset_gisaid_subcol.loc[:,'rss_exposure'] = subset_gisaid_subcol.country_exposure
subset_gisaid_subcol.loc[:,'rta_exposure'] = subset_gisaid_subcol.country_exposure
#print(subset_gisaid_subcol)

#print(subset_lspq_subcol.columns)
#print(subset_gisaid_subcol.columns)

subset_lspq_subcol.to_csv(out_metadata,sep="\t",index=False)
subset_gisaid_subcol.to_csv(out_metadata,sep="\t",index=False,mode='a',header=False)

SeqIO.write(rec_list,out_sequences,'fasta')

logging.info("End")
