"""
Eric Fournier 2020-04-14

"""

import pypostalcode
import pycountry

LEVEL = {'first':'rta','second':'province','third':'country'}
 
pcdb = pypostalcode.PostalCodeDatabase()

def GetCanadaRTA():
    CANADA_RTA = []

    for loc_obj in pcdb.find_postalcode(province = None):
        CANADA_RTA.append(loc_obj.postalcode)

    return CANADA_RTA

def GetQuebecRTA():
    QUEBEC_RTA = []
    
    for loc_obj in pcdb.find_postalcode(province = 'Quebec'):
        QUEBEC_RTA.append(loc_obj.postalcode)

    return QUEBEC_RTA


def GetCanadaProvinces():
    pypostalcode_conn = pypostalcode.ConnectionManager()
    sql = "SELECT distinct province from PostalCodes"

    CANADA_PROV = [p[0] for p in pypostalcode_conn.query(sql)]

    return CANADA_PROV

def GetCountries():
    COUNTRIES = []
    
    for c in pycountry.countries:
        COUNTRIES.append(c.name)

    return COUNTRIES
