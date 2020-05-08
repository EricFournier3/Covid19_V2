"""
Eric Fournier 2020-04-22
"""

import os

base_dir =  "/data/PROJETS/COVID_19/TestMetadata/"
virus = "ncov"


file_in = os.path.join(base_dir,"export.tsv")
file_in_handle = open(file_in)
file_in_handle.readline()

metadata = os.path.join(base_dir,"metadata.tsv")
metadata_handle = open(metadata,'w')

metadata_handle.write("strain\t" + "virus\t" + "date\t" +  "age\t" + "sex\t" + "rss\t" + "rta\t" + "rss_exposure\t" + "rta_exposure\n")

for line in file_in_handle:
    s = (line.strip('\n'))
    mylist = s.split('\t')
    #print(str(mylist))
    no_lspq = mylist[0]
    age = mylist[1]
    sex = mylist[2]
    rss = mylist[3]
    voyage = mylist[5]
    date_prelev = mylist[6]
    rta = mylist[9][0:3]

    if voyage == "AUCUN_VOYAGE":
        rss_exposure = rss
        rta_exposure = rta
    else:
        rss_exposure = voyage
        rta_exposure = voyage

    metadata_handle.write(str(no_lspq) + "\t" + virus + "\t" + str(date_prelev) + "\t" + str(age) + "\t" + str(sex) +
    "\t" + str(rss) + "\t" + str(rta) + "\t"  + str(rss_exposure) +  "\t" + str(rta_exposure) + "\n" )

file_in_handle.close()

metadata_handle.close()
