import argparse
import pandas as pd
import yaml
import os
import subprocess


# Forced colours MUST NOT appear in the ordering TSV
forced_colors = {
}

def CleanUserColorsDict(TraitName,user_color):

    cleaned_user_color = {}

    rta_lat_long_file = os.path.join(os.getcwd(),"config/rta_lat_long.tsv")
    rss_lat_long_file = os.path.join(os.getcwd(),"config/rss_lat_long.tsv")
   
    lookup = {'rss':rss_lat_long_file,'rta':rta_lat_long_file}
 
    if (TraitName in ['country','division']):
        for my_trait in user_color.keys():

            check_rta = 0
            check_rss = 0
            try:
                check_rta = subprocess.check_output(['grep','-c',my_trait, lookup['rta']])
            except:
                pass
           
            try:
                check_rss = subprocess.check_output(['grep','-c',my_trait, lookup['rss']])
            except:
                pass

            
            if (int(check_rta) > 0) or (int(check_rss) > 0):
                continue
            else:
                cleaned_user_color[my_trait] = user_color[my_trait] 
       
    elif (TraitName == 'rss'): 
        for my_trait in user_color.keys():

            check_rta = 0

            try:
                check_rta = subprocess.check_output(['grep','-c',my_trait, lookup['rta']])
            except:
                pass
            
            if (int(check_rta) > 0):
                continue
            else:
                cleaned_user_color[my_trait] = user_color[my_trait] 

    elif (TraitName == 'rta'): 
        for my_trait in user_color.keys():
            check_rss = 0

            try:
                check_rss = subprocess.check_output(['grep','-c',my_trait, lookup['rss']])
            except:
                pass

            if (int(check_rss) > 0):
                continue
            else:
                cleaned_user_color[my_trait] = user_color[my_trait] 
        
    return cleaned_user_color


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Assign colors based on ordering",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--ordering', type=str, required=True, help="input ordering file")
    parser.add_argument('--color-schemes', type=str, required=True, help="input color schemes file")
    parser.add_argument('--metadata', type=str, help="if provided, restrict colors to only those found in metadata")
    parser.add_argument('--output', type=str, required=True, help="output colors tsv")

    parser.add_argument('--user-colors', type=str, required=True, help="yaml input for user defined color")
    args = parser.parse_args()

    user_colors = args.user_colors
    focal = True
    color_handle = open(user_colors)
    color_yaml = yaml.load(color_handle)
    user_colors_dict = color_yaml['colors'][0] 

    assignment = {}
    with open(args.ordering) as f:
        for line in f.readlines():
            array = line.lstrip().rstrip().split("\t")
            if len(array) == 2:
                name = array[0]
                trait = array[1]
                if name not in assignment:
                    assignment[name] = [trait]
                else:
                    assignment[name].append(trait)

    # if metadata supplied, go through and
    # 1. remove assignments that don't exist in metadata
    # 2. remove assignments that have 'focal' set to 'False' in metadata
    if args.metadata:
        metadata = pd.read_csv(args.metadata, delimiter='\t')
        for name, trait in assignment.items():
            if name in metadata:
                subset_present = [x for x in assignment[name] if x in metadata[name].unique()]
                assignment[name] = subset_present
            if name in metadata and 'focal' in metadata:
                focal_list = metadata.loc[metadata['focal'] == True, name].unique()
                subset_focal = [x for x in assignment[name] if x in focal_list]
                assignment[name] = subset_focal
                focal = True

    schemes = {}
    counter = 0
    with open(args.color_schemes) as f:
        for line in f.readlines():
            counter += 1
            array = line.lstrip().rstrip().split("\t")
            schemes[counter] = array

    with open(args.output, 'w') as f:
        for trait_name, trait_array in assignment.items():

            if focal:
                cleaned_user_colors_dict = CleanUserColorsDict(trait_name,user_colors_dict)
                extra_trait_values = list(cleaned_user_colors_dict.keys())
                extra_color_values = list(cleaned_user_colors_dict.values())

            else:
                extra_trait_values = []
                extra_color_values = []

            #trait_array = set(trait_array)

            extra_trait_values_set = set(extra_trait_values)

            trait_array = list(set(trait_array) - extra_trait_values_set)

            if len(trait_array)==0 and len(extra_trait_values)==0:
                print(f"No traits found for {trait_name}")
                continue
            try:
                color_array = schemes[len(trait_array)-1]
            except:
                color_array = []

            zipped = list(zip(trait_array+extra_trait_values, color_array+extra_color_values))
            for trait_value, color in zipped:
                f.write(trait_name + "\t" + trait_value + "\t" + color + "\n")
            f.write("\n")
