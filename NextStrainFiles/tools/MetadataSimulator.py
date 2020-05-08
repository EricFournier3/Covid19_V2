"""
Eric Fournier 20200415
"""

import os
from SequenceSimulator import CreateFastaCovid
import random
import datetime


location_file = os.path.join(os.path.dirname(os.getcwd()),"config/ordering.tsv")

class RandomDataGenerator(object):
    def __init__(self):
        pass

    def CreateRandomList(self, target_length, data_list):
        rnd_list = [] 

        for i in range(0,target_length):
            rnd_index = random.randint(0,len(data_list) - 1)
            rnd_list.append(data_list[rnd_index]) 

        return(rnd_list.copy())

    def CreateRandomDateList(self,target_length):
        start_date = datetime.date(year=2018, month=1, day=1)
        end_date = datetime.date(year=2020, month=11, day=1)
        rnd_date_list = []

        for j in range(0, target_length):
            time_between_dates = end_date - start_date
            days_between_dates = time_between_dates.days
            random_number_of_days = random.randrange(days_between_dates)
            random_date = start_date + datetime.timedelta(days=random_number_of_days)
            rnd_date_list.append(str(random_date))

        return(rnd_date_list.copy())

class MetadataBuilder(object):
    def __init__(self, number_of_isolate, is_homogeneous, fasta_in, fasta_out, metadata):
        self.fasta_in = fasta_in
        self.fasta_out = fasta_out

        self.location_file = os.path.join(os.path.dirname(os.getcwd()),"config/ordering.tsv")
        self.number_of_isolate = number_of_isolate
        self.is_homogegneous = is_homogeneous
        self.metadata = metadata

        self.header = {'Strain' : 'strain', 'Virus' : 'virus', 'Date' : 'date' ,'Rta' : 'rta', 'CountryExposure' : 'country_exposure',
        'ProvinceExposure' : 'province_exposure', 'Rta_exposure' : 'rta_exposure', 'Segment' : 'segment', 'Length' : 'length',
        'Host' : 'host', 'Age' : 'age', 'Sex' : 'sex', 'LabOrig' : 'originating_lab', 'LabSubmt' : 'submitting_lab', 'Authors' : 'authors',
        'Url' : 'url', 'Title' : 'title'}

        self.static_data = {'virus' : 'ncov', 'segment' : 'genome', 'host' : 'human', 'age' : '30',  'sex' : 'Male', 'originating_lab' : 'Labx',
        'submitting_lab' : 'Laby', 'authors' : 'Somebody', 'url' : 'NA', 'title' : 'NA'}

        self.location_dict = {'rta' : set(), 'province' : set(), 'country' : set()}
        self.rta_list = []
        self.province_list = []
        self.country_list = []

        self.rnd_generator = RandomDataGenerator()

        self.random_rta_list = None
        self.random_rta_exposure_list = None
        self.random_province_list = None
        self.random_country_list = None
        self.random_date_list = None

        self.strain_list,self.seq_len_list = CreateFastaCovid(fasta_in_file = self.fasta_in, fasta_out_file = self.fasta_out, nb_seq = self.number_of_isolate, homogeneous = is_homogeneous)
        self.BuildLocationList()
        self.CreateRandomData()

    def CreateRandomData(self):
        self.random_rta_list = self.rnd_generator.CreateRandomList(self.number_of_isolate,self.rta_list)
        self.random_rta_exposure_list = self.rnd_generator.CreateRandomList(self.number_of_isolate,self.rta_list)
        self.random_province_list = self.rnd_generator.CreateRandomList(self.number_of_isolate,self.province_list)
        self.random_country_list = self.rnd_generator.CreateRandomList(self.number_of_isolate,self.country_list)
        self.random_date_list = self.rnd_generator.CreateRandomDateList(self.number_of_isolate)

    def BuildLocationList(self): 
        with open(self.location_file) as readf:
            for line in readf:
                info = line.split('\t')
                self.location_dict[info[0]].add(info[1].strip('\n'))

        self.rta_list = list(self.location_dict['rta'] - self.location_dict['province'].union(self.location_dict['country']))
        self.province_list = list(self.location_dict['province'])
        self.country_list = list(self.location_dict['country'])

    def CreateMetadataFile(self):
        with open(self.metadata,'w') as metadata_handle:
            metadata_handle.write(self.header['Strain'] + "\t" + self.header['Virus'] + "\t" + self.header['Date'] + "\t" + self.header['Rta'] +
            "\t" + self.header['CountryExposure'] + "\t" + self.header['ProvinceExposure'] + "\t" + self.header['Rta_exposure'] +
            "\t" + self.header['Segment'] + "\t" + self.header['Length'] + "\t" + self.header['Host'] + "\t" + self.header['Age'] +
            "\t" + self.header['Sex'] + "\t" + self.header['LabOrig'] + "\t" + self.header['LabSubmt'] + "\t" + self.header['Authors'] +
            "\t" + self.header['Url'] + "\t" + self.header['Title'] + "\n")

            index = 0

            for strain in self.strain_list:
                metadata_handle.write(strain + "\t")
                metadata_handle.write(self.static_data['virus'] + "\t")
                metadata_handle.write(self.random_date_list[index] + "\t")
                metadata_handle.write(self.random_rta_list[index] + "\t")
                metadata_handle.write(self.random_country_list[index] + "\t")
                metadata_handle.write(self.random_province_list[index] + "\t")
                metadata_handle.write(self.random_rta_exposure_list[index] + "\t")
                metadata_handle.write(self.static_data['segment'] + "\t")
                metadata_handle.write(self.seq_len_list[index] + "\t")
                metadata_handle.write(self.static_data['host'] + "\t")
                metadata_handle.write(self.static_data['age'] + "\t")
                metadata_handle.write(self.static_data['sex'] + "\t")
                metadata_handle.write(self.static_data['originating_lab'] + "\t")
                metadata_handle.write(self.static_data['submitting_lab'] + "\t")
                metadata_handle.write(self.static_data['authors'] + "\t")
                metadata_handle.write(self.static_data['url'] + "\t")
                metadata_handle.write(self.static_data['title'])
                metadata_handle.write('\n')

                index += 1


fasta_in_file =  os.path.join(os.path.dirname(os.getcwd()),"data/gisaid_covid19_usa_1759.fasta")

fasta_out =  os.path.join(os.path.dirname(os.getcwd()),"data/covid19_10seq_homogeneous.fasta")
metadata_file = os.path.join(os.path.dirname(os.getcwd()),"data/metadata_10seq_homogeneous.tsv")

number_of_isolate = 10
is_homogeneous_seq = True

mbuilder = MetadataBuilder(number_of_isolate, is_homogeneous_seq, fasta_in_file, fasta_out, metadata_file)
mbuilder.CreateMetadataFile()
