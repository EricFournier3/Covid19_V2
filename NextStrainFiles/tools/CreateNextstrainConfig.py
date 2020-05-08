"""
Eric Fournier 2020-04-14

"""


import Locations
import pickle
import logging
import ssl
import pgeocode
import os
from geopy.geocoders import Nominatim

#see https://stackoverflow.com/questions/50236117/scraping-ssl-certificate-verify-failed-error-for-http-en-wikipedia-org
ssl._create_default_https_context = ssl._create_unverified_context
nomi = pgeocode.Nominatim('ca')

geolocator = Nominatim(user_agent = "my_application")
logging.basicConfig(level = logging.DEBUG)

#CANADA_RTA = Locations.GetCanadaRTA()
QUEBEC_RTA = Locations.GetQuebecRTA()
CANADA_PROV = Locations.GetCanadaProvinces()
COUNTRIES = Locations.GetCountries()

QUEBEC_LAT_LONG = {}
CANADA_PROV_LAT_LONG = {}
COUNTRIES_LAT_LONG = {}

quebec_lat_long_obj = os.path.join(os.path.dirname(os.getcwd()),"LatLongObj/quebec_lat_long.pkl")
canada_prov_lat_long_obj = os.path.join(os.path.dirname(os.getcwd()),"LatLongObj/canada_prov_lat_long.pkl")
countries_lat_long_obj = os.path.join(os.path.dirname(os.getcwd()),"LatLongObj/countries_lat_long.pkl")

lat_long_out = os.path.join(os.path.dirname(os.getcwd()),"config/lat_longs.tsv")
ordering_out = os.path.join(os.path.dirname(os.getcwd()),"config/ordering.tsv")

level_1 = Locations.LEVEL['first']
level_2 = Locations.LEVEL['second']
level_3 = Locations.LEVEL['third']

def GetQuebecRTA_LatLong():

    logging.info(">>>> Get Quebec RTA Latitude Longitude")

    for rta in QUEBEC_RTA:
        try:
            #logging.info("Try for " + rta)
            latitude = Locations.pcdb[rta].latitude
            longitude = Locations.pcdb[rta].longitude
            temp_dict = {'latitude': str(latitude), 'longitude':str(longitude)}
            QUEBEC_LAT_LONG[rta] = temp_dict.copy() 
            #print(QUEBEC_LAT_LONG[rta])
        except Exception as e:
            logging.error("**************      Error with rta " + rta + " >> " + e.__str__())

        out = open(quebec_lat_long_obj,"wb")
        pickle.dump(QUEBEC_LAT_LONG,out,-1)
        out.close()

def GetCanadaProvinces_LatLong():
    
    logging.info(">>>> Get Canada Provinces Latitude Longitude")
    
    for prov in CANADA_PROV:
        try:
            loc = geolocator.geocode(prov)
            latitude = loc.latitude
            longitude = loc.longitude 
            temp_dict = {'latitude': str(latitude), 'longitude': str(longitude)}
            CANADA_PROV_LAT_LONG[prov] = temp_dict.copy()
        except :
            logging.error("***************   Error with province " + prov)

        out = open(canada_prov_lat_long_obj,"wb")
        pickle.dump(CANADA_PROV_LAT_LONG,out,-1)
        out.close()

def GetCountries_LatLong():

    logging.info(">>>> Get Countries Latitude Longitude")
     
    for country in COUNTRIES:
        try:
            loc = geolocator.geocode(country)
            latitude = loc.latitude
            longitude = loc.longitude 
            temp_dict = {'latitude': str(latitude), 'longitude': str(longitude)}
            COUNTRIES_LAT_LONG[country] = temp_dict.copy()
        except :
            logging.error("**************    Error with country " + country)

        out = open(countries_lat_long_obj,"wb")
        pickle.dump(COUNTRIES_LAT_LONG,out,-1)
        out.close()

if not os.path.exists(quebec_lat_long_obj):
    GetQuebecRTA_LatLong()
else:
    QUEBEC_LAT_LONG = pickle.load(open(quebec_lat_long_obj,"rb"))

if not os.path.exists(canada_prov_lat_long_obj):
    GetCanadaProvinces_LatLong()
else:
    CANADA_PROV_LAT_LONG = pickle.load(open(canada_prov_lat_long_obj,"rb"))

if not os.path.exists(countries_lat_long_obj):
    GetCountries_LatLong()
else:
    COUNTRIES_LAT_LONG = pickle.load(open(countries_lat_long_obj, "rb"))


def WriteQuebecRTA_LatLong(handler):
    
    for rta in QUEBEC_LAT_LONG.keys():
        handler.write(level_1 + "\t" + rta + "\t" + str(format(float(QUEBEC_LAT_LONG[rta]['latitude']),'.6f')) + "\t" + str(format(float(QUEBEC_LAT_LONG[rta]['longitude']),'.6f')) + "\n" )

def WriteCanadaProvince_LatLong(handler):

    for prov in CANADA_PROV_LAT_LONG.keys():
        handler.write(level_1 + "\t" + prov + "\t" + str(format(float(CANADA_PROV_LAT_LONG[prov]['latitude']),'.6f')) + "\t" + str(format(float(CANADA_PROV_LAT_LONG[prov]['longitude']),'.6f')) + "\n" )

    for prov in CANADA_PROV_LAT_LONG.keys():
        handler.write(level_2 + "\t" + prov + "\t" + str(format(float(CANADA_PROV_LAT_LONG[prov]['latitude']),'.6f')) + "\t" + str(format(float(CANADA_PROV_LAT_LONG[prov]['longitude']),'.6f')) + "\n" )


def WriteContries_LatLong(handler):

    for country in COUNTRIES_LAT_LONG.keys():
        handler.write(level_1 + "\t" + country + "\t" + str(format(float(COUNTRIES_LAT_LONG[country]['latitude']),'.6f')) + "\t" + str(format(float(COUNTRIES_LAT_LONG[country]['longitude']),'.6f')) + "\n" )

    for country in COUNTRIES_LAT_LONG.keys():
        handler.write(level_3 + "\t" + country + "\t" + str(format(float(COUNTRIES_LAT_LONG[country]['latitude']),'.6f')) + "\t" + str(format(float(COUNTRIES_LAT_LONG[country]['longitude']),'.6f')) + "\n" )

def CreateLatitudesLongitudesFile():
    writer = open(lat_long_out,'w') 

    WriteQuebecRTA_LatLong(handler = writer)
    WriteCanadaProvince_LatLong(handler = writer)
    WriteContries_LatLong(handler = writer)

    writer.close()

def WriteRTA(handler):
    for rta in QUEBEC_LAT_LONG.keys():
        handler.write(level_1 + "\t" + rta + "\n")

    for prov in CANADA_PROV_LAT_LONG.keys():
        handler.write(level_1 + "\t" + prov + "\n")

    for country in COUNTRIES_LAT_LONG.keys():
        handler.write(level_1 + "\t" + country + "\n")

def WriteProvince(handler):
    for prov in CANADA_PROV_LAT_LONG.keys():
        handler.write(level_2 + "\t" + prov + "\n")

def WriteCountries(handler):
    for country in COUNTRIES_LAT_LONG.keys():
        handler.write(level_3 + "\t" + country + "\n")

def CreateOrderingFile():
    writer = open(ordering_out,'w')
    
    WriteRTA(writer)
    WriteProvince(writer)
    WriteCountries(writer)

    writer.close()

CreateLatitudesLongitudesFile()
CreateOrderingFile()
