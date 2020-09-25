# coding=utf-8
"""
Eric Fournier 2020-09-24

"""


import mysql.connector
import pandas as pd

class MySQLcovid19:
    host = 'localhost'
    user = 'root'
    password = 'lspq2019'
    database = 'TestCovid19_20200923'
    connection = None

    @classmethod
    def SetConnection(cls):
        cls.connection =  mysql.connector.connect(host=cls.host,user=cls.user,password=cls.password,database=cls.database)
  
    @classmethod
    def GetConnection(cls):
        return cls.connection

    @classmethod
    def GetCursor(cls):
        return cls.GetConnection().cursor()

    @classmethod
    def Commit(cls):
        cls.connection.commit()


class MySQLcovid19Selector:

    @classmethod
    def GetSampleDate(cls,cursor,sample_name):

        cursor.execute("select DATE_PRELEV_HOPITAL from Prelevements where GENOME_QUEBEC_REQUETE  like '%{0}%'".format(sample_name))
        #date_prelev = cursor.fetchone()
        date_prelev = cursor.fetchall()

        if len(date_prelev) == 1 :
            date_prelev = date_prelev[0][0]
            mois = Utils.GetFrenchMonth(date_prelev.strftime("%B"))
            return(date_prelev.strftime("%d {0} %Y".format(mois)))

        return("indéterminé")

    @classmethod
    def GetMetadataAsPdDataFrame(cls,conn,spec_list,metadata_destination,extract_all_samples):
        import time
        
        start = time.time()

        columns_renamed = {'GENOME_QUEBEC_REQUETE':'sample','DATE_PRELEV_HOPITAL':'sample_date','TRAVEL_HISTORY':'country_exposure','CT':'ct','RSS_LSPQ_CAS':'rss','SEXEINFO':'sex','COUNTRY':'country','DIVISION':'division','DTNAISSINFO':'date_naiss'}        

        spec_list = '|'.join(spec_list)
        prelevements_alias = 'pr'
        patients_alias = 'p'

        PRELEVEMENTS_COLUMNS = ['GENOME_QUEBEC_REQUETE','DATE_PRELEV_HOPITAL','TRAVEL_HISTORY','CT']
        PRELEVEMENTS_COLUMNS = [prelevements_alias + "." + col for col in PRELEVEMENTS_COLUMNS]
        PRELEVEMENTS_COLUMNS = ','.join(PRELEVEMENTS_COLUMNS)

        if (metadata_destination == 'lspq'):
            PATIENTS_COLUMNS = ['RSS_LSPQ_CAS','DTNAISSINFO','SEXEINFO']
        else:
            PATIENTS_COLUMNS = ['SEXEINFO']

        PATIENTS_COLUMNS = [patients_alias + "." + col for col in PATIENTS_COLUMNS]
        PATIENTS_COLUMNS = ','.join(PATIENTS_COLUMNS)
        COUNTRY = "'Canada'"
        DIVISION = "'Quebec'"
        RTA = "'G8P'"
        if(extract_all_samples):
            sql = "SELECT {0},{1}, {3} as COUNTRY, {4} as DIVISION FROM Prelevements pr inner join Patients p on  p.ID_PATIENT = pr.ID_PATIENT ".format(PRELEVEMENTS_COLUMNS,PATIENTS_COLUMNS,spec_list,COUNTRY,DIVISION)
        else:
            sql = "SELECT {0},{1}, {3} as COUNTRY, {4} as DIVISION FROM Prelevements pr inner join Patients p on  p.ID_PATIENT = pr.ID_PATIENT WHERE pr.GENOME_QUEBEC_REQUETE REGEXP '{2}' ".format(PRELEVEMENTS_COLUMNS,PATIENTS_COLUMNS,spec_list,COUNTRY,DIVISION)

        df = pd.read_sql(sql,con=conn)
        df = df.rename(columns=columns_renamed)

        end = time.time()

        print("In GetMetadataAsPdDataFrame for ",end - start, " seconds")
        return df


class Utils():

    @classmethod
    def GetFrenchMonth(cls,english_month):
        month_map = {'January':'Janvier','February':'Février','March':'Mars','April':'Avril','May':'Mai','June':'Juin','Juillet':'July','August':'Août',
                     'September':'Septembre','October':'Octobre','November':'Novembre','December':'Décembre'} 

        return(month_map[english_month])
