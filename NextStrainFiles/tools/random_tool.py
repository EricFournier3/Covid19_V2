import pandas
import os
import pickle

#TODO Pandas
in_dat_base_dir = "/data/Applications/GitScript/Covid19/NextStrainFiles/"

rta_lat_long_obj = os.path.join(in_dat_base_dir,"python_obj","quebec_lat_long.pkl")
province_lat_long_obj = os.path.join(in_dat_base_dir,"python_obj","canada_prov_lat_long.pkl")

rss_lat_long = {"01-Bas-Saint-Laurent" : {"latitude" : 48.655230,"longitude" : -67.467500}, "02-Saguenay - Lac-Saint-Jean" : {"latitude":49.378300,"longitude":-72.730600},"03-Capitale-Nationale" : {"latitude":47.376000,"longitude":-71.123370}, "04-Mauricie et Centre-du-Québec" : {"latitude":46.663000,"longitude":-72.851200},
                "04-Mauricie et Centre-du-Québec" : {"latitude":46.663000,"longitude" : -72.851200}, "05-Estrie" : {"latitude":45.404200,"longitude":-71.892900}, "06-Montréal" : {"latitude" : 45.501700,"longitude" : -73.567300}, "07-Outaouais" : {"latitude":45.476500,"longitude":-75.701300}, "08-Abitibi-Témiscamingue" : {"latitude":47.734400,"longitude":-77.672800}, "09-Côte-Nord" : {"latitude":50.8967,"longitude":-66.149700},"10-Nord-du-Québec" : {"latitude":52.739900,"longitude": -73.549100}, "11-Gaspésie - Îles-de-la-Madeleine" : {"latitude":47.573100, "longitude":-62.325500}, "12-Chaudière-Appalaches" : {"latitude":46.439000,"longitude":-71.021300}, "13-Laval" : {"latitude":45.606600,"longitude": -73.712400}, "14-Lanaudière" : {"latitude":46.125600,"longitude" : -73.704200}, "15-Laurentides" : {"latitude":46.618200,"longitude":-75.018100}, "16-Montérégie" : {"latitude":45.501700,"longitude" : -73.567300}, "17-Nunavik" : {"latitude" : 60.034200,"longitude" : -70.011800}, "18-Terres-Cries-de-la-Baie-James" : {"latitude":51.700800,"longitude" : -76.046800}}


lat_long_gisaid = os.path.join(in_dat_base_dir,"config/original/config/lat_longs.tsv")


LEVEL = {"first" : "country", "second" : "division" , "third" : "rss", "fourth" : "rta"}

#print(str(pickle.load(open(province_lat_long_obj,'rb'))))

with open(os.path.join(in_dat_base_dir,"config","province_lat_long.tsv"),'w') as wr:
    for key, val in pickle.load(open(province_lat_long_obj,'rb')).items():
        print("KEY ",key, str(val))
        wr.write(key + "\t" + str(format(float(val["latitude"]),'.6f')    ) + "\t" + str(format(float(val["longitude"]),'.6f')) + "\n")
        
with open(os.path.join(in_dat_base_dir,"config","rta_lat_long.tsv"),'w') as wr:
    for key, val in pickle.load(open(rta_lat_long_obj,'rb')).items():
        print("KEY ",key, str(val))
        wr.write(key + "\t" + str(format(float(val["latitude"]),'.6f')    ) + "\t" + str(format(float(val["longitude"]),'.6f')) + "\n")

with open(os.path.join(in_dat_base_dir,"config","rss_lat_long.tsv"),'w') as wr:
    for key, val in rss_lat_long.items():
        print("KEY ",key, str(val))
        wr.write(key + "\t" + str(format(float(val["latitude"]),'.6f')    ) + "\t" + str(format(float(val["longitude"]),'.6f')) + "\n")
