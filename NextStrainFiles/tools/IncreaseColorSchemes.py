import os
import shutil


small_color_schemes = os.path.join(os.path.dirname(os.getcwd()),"config/color_schemes_small.tsv")
large_color_schemes = os.path.join(os.path.dirname(os.getcwd()),"config/color_schemes.tsv")

shutil.copyfile(small_color_schemes,large_color_schemes)

schemes_list = open(small_color_schemes).readlines()
nb_schemes = len(schemes_list)
max_schemes = schemes_list[nb_schemes - 1].split('\t')

large_color_schemes_handle = open(large_color_schemes,'a')
small_color_schemes_handle = open(small_color_schemes)

for line in small_color_schemes_handle:
    arr = line.split('\t')
    new_schemes = max_schemes + arr
    new_schemes = '\t'.join(new_schemes)
    large_color_schemes_handle.write(str(new_schemes))


small_color_schemes_handle.close()
large_color_schemes_handle.close()
