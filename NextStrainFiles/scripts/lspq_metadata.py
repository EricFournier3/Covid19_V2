# -*- coding: utf-8 -*-


import pandas as pd
import datetime

def from_dob_to_age(born):
    today = datetime.date.today()
    return today.year - born.year - ((today.month, today.day) < (born.month, born.day))


def GetQcMetadataDataFrame(qc_seq_id,metadata_lspq):
    renamed_columns_dict = {'sample':'strain','sample_date':'date','date_naiss':'DTNAISSINFO'}
    #pd_qc_metadata = pd.read_csv("/home/foueri01@inspq.qc.ca/temp/TEST_GET_FASTA_FOR_NEXTSTRAIN/OUT/METADATA/metadata_2020-10-21_from_updateListOfRuns_v2_2020-10-21_consensusList.list.tsv",sep="\t",index_col=False,dtype={'ct':str,'rta':str})
    pd_qc_metadata = pd.read_csv(metadata_lspq,sep="\t",index_col=False,dtype={'ct':str,'rta':str})

    #Modif_20210126 on desactive
    #pd_qc_metadata = pd_qc_metadata.loc[pd_qc_metadata['sample'].isin(qc_seq_id),:]

    #Modif_20210126 add this line
    pd_qc_metadata['sample'] = pd_qc_metadata['fasta_id']

    #print("************************** SHAPE 1 ",pd_qc_metadata.shape[0])

    pd_qc_metadata = pd_qc_metadata.rename(columns=renamed_columns_dict)
    pd_qc_metadata['virus'] = 'ncov'


    pd_qc_metadata.loc[pd_qc_metadata['country_exposure'].isin(['AUCUN_VOYAGE','INDETERMINE']),'rta_exposure'] = pd_qc_metadata['rta']
    #print(pd_qc_metadata)
    pd_qc_metadata.loc[pd_qc_metadata['country_exposure'].isin(['AUCUN_VOYAGE','INDETERMINE']),'rss_exposure'] = pd_qc_metadata['rss']
    pd_qc_metadata.loc[pd_qc_metadata['country_exposure'].isin(['AUCUN_VOYAGE','INDETERMINE']),'division_exposure'] = pd_qc_metadata['division']
    
    pd_qc_metadata.loc[pd_qc_metadata['country_exposure'].str.upper().isin(['CANADA']),'country_exposure'] = 'Canada'
   
 
    pd_qc_metadata.loc[~pd_qc_metadata['country_exposure'].isin(['AUCUN_VOYAGE','INDETERMINE']),'rta_exposure'] = pd_qc_metadata['country_exposure']
    pd_qc_metadata.loc[~pd_qc_metadata['country_exposure'].isin(['AUCUN_VOYAGE','INDETERMINE']),'rss_exposure'] = pd_qc_metadata['country_exposure']
    pd_qc_metadata.loc[~pd_qc_metadata['country_exposure'].isin(['AUCUN_VOYAGE','INDETERMINE']),'division_exposure'] = pd_qc_metadata['country_exposure']
    
    pd_qc_metadata.loc[pd_qc_metadata['country_exposure'].isin(['AUCUN_VOYAGE','INDETERMINE']),'country_exposure'] = pd_qc_metadata['country']

    now = pd.to_datetime('today')
    pd_qc_metadata['DTNAISSINFO'] = pd.to_datetime(pd_qc_metadata['DTNAISSINFO'])
    pd_qc_metadata['age'] = pd_qc_metadata['DTNAISSINFO'].apply(lambda x: from_dob_to_age(x))

    #TODO CE NE SERA PAS TOUJOURS 2020 => BESOIN DE LA MEME ANNEE QUE CELLE QUI EST DANS HEADER FASTA
    #pd_qc_metadata["strain"] = "Canada/Qc-" + pd_qc_metadata['strain'] + "/" + pd_qc_metadata['date'].str[:4]

    #Modif_20210126 on desactive
    #pd_qc_metadata["strain"] = "Canada/Qc-" + pd_qc_metadata['strain'] + "/2020"

    pd_qc_metadata['rss'] =  pd_qc_metadata['rss'].str.replace('é','e').str.replace('ô','o').str.replace('î','i').str.replace('è','e').str.replace('Î','i').str.replace(' - ','-')
    pd_qc_metadata['rss_exposure'] =  pd_qc_metadata['rss_exposure'].str.replace('é','e').str.replace('ô','o').str.replace('î','i').str.replace('è','e').str.replace('Î','i').str.replace(' - ','-')

    pd_qc_metadata = pd_qc_metadata[['strain','virus','date','sex','age','country','division','rss','rta','country_exposure','division_exposure','rss_exposure','rta_exposure','ct','OUTBREAK']]
    pd_qc_metadata.loc[pd_qc_metadata['ct'] == '0','ct'] = '99'
    pd_qc_metadata.loc[pd_qc_metadata['OUTBREAK'].isnull(),['OUTBREAK']] = 'NoOutbreakRelated'
    #print("************************** SHAPE 2 ",pd_qc_metadata.shape[0])
    mystrain = pd_qc_metadata['strain']
    test_df = pd_qc_metadata[mystrain.isin(mystrain[mystrain.duplicated()])]
    #print(test_df)

    pd_qc_metadata = pd_qc_metadata.drop_duplicates(subset='strain',keep='first')
    #print("************************** SHAPE 3 ",pd_qc_metadata.shape[0])
    #print(pd_qc_metadata)
    return(pd_qc_metadata)
