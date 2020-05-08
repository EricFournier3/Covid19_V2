import os, json, sys
import pypostalcode
import pycountry

'''
json_file = "/data/PROJETS/COVID_19/Test4/config/auspice_config.json"

try:
    with open(json_file) as ifile:
        config = json.load(ifile)
except json.decoder.JSONDecodeError as err:
    print("\tError message: '{}'".format(err.msg))
    print("\tLine number: '{}'".format(err.lineno))
    print("\tColumn number: '{}'".format(err.colno))
'''

pcdb = pypostalcode.PostalCodeDatabase()
#print(pcdb['G1Y'].latitude)
i = 0
for loc_obj in pcdb.find_postalcode(province = 'Quebec'):
    #print(loc_obj.postalcode)
    i+=1
print(i)
