#voir https://www.latlong.net/category/provinces-40-60.html

import pgeocode
import pypostalcode
import pypostalcode
import ssl
import logging
import geopy

logging.basicConfig(level=logging.DEBUG)
#voir https://stackoverflow.com/questions/50236117/scraping-ssl-certificate-verify-failed-error-for-http-en-wikipedia-org
ssl._create_default_https_context = ssl._create_unverified_context
nomi = pgeocode.Nominatim('ca')
#print(nomi.query_postal_code('J7E').latitude)

prov_lat_long = {'Alberta':{'lat':'55.000000','long':'-115.000000'},'Manitoba':{'lat':'53.760860','long':'-98.813873'},
'New Brunswick':{'lat':'46.498390','long':'-66.159668'},'Newfoundland and Labrador':{'lat':'53.135509','long':'-57.660435'},
'Newfoundland and Labrador':{'lat':'53.135509','long':'-57.660435'},'Northwest Territory':{'lat':'64.825500','long':'-124.845700'},
'Nova Scotia':{'lat':'45.000000','long':'-63.000000'},'Nunavut Territory':{'lat':'70.453262','long':'-86.798981'},
'Ontario':{'lat':'50.000000','long':'-85.000000'},'Prince Edward Island':{'lat':'46.250000','long':'-63.000000'},'Quebec':{'lat':'53.000000','long':'-70.000000'},'Saskatchewan':{'lat':'55.000000','long':'-106.000000'},'Yukon':{'lat':'64.282300','long':'135.000000'},'British Columbia':{'lat':'53.726669','long':'-127.647621'}}

pcdb = pypostalcode.PostalCodeDatabase()

loc_canada = []
loc_quebec = []

#for loc in pcdb.find_postalcode(province=None):
#    loc_canada.append(loc)

for loc in pcdb.find_postalcode(province='Quebec'):
    loc_quebec.append(loc)
    #print(loc.postalcode)

myloc = pcdb['J7E']
#print(myloc.latitude, ' ',myloc.longitude)

#print('Canada ' + str(len(loc_canada)))
#print('Quebec ' + str(len(loc_quebec)))

lat_long_file = '/data/Applications/GitScript/Covid19/LatLong.txt'

lat_long_file_handler = open(lat_long_file,'w')
'''

for loc in loc_quebec:
    rta = loc.postalcode
    logging.info('Write latitude longitude for ' + rta)
    lat_long_file_handler.write(rta + '\t' + str(format(pcdb[rta].latitude,'.6f')) + '\t' + str(format(pcdb[rta].longitude,'.6f')) + '\n')
'''



for prov in prov_lat_long.keys():
    logging.info("Write latitude longitude for " + prov)
    lat_long_file_handler.write(prov + "\t" + prov_lat_long[prov]['lat'] + "\t" + prov_lat_long[prov]['long'] + '\n')


lat_long_file_handler.close()

logging.info('Finish')

myconn = pypostalcode.ConnectionManager()
PC_FIND_ALL_PROVINCE = "SELECT distinct province from PostalCodes"
prov_tuple = myconn.query(PC_FIND_ALL_PROVINCE)
#print(prov_tuple)
prov_list = [x[0] for x in prov_tuple]
#print(prov_list)

from geopy.geocoders import Nominatim
geolocator = Nominatim(user_agent="my_application")

lat_long_file_handler = open(lat_long_file,'w')

for prov in prov_list:
    print("PROV IS " + prov)
    try:
        location = geolocator.geocode(prov)
        logging.info("Write latitude longitude for " + prov)
        lat_long_file_handler.write(prov + "\t" + str(location.latitude) + "\t" + str(location.longitude) + '\n')
    except:
        logging.error("No lat/long for " + prov)

lat_long_file_handler.close()
