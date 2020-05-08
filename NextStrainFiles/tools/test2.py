import os, json, sys

json_file = "/data/PROJETS/COVID_19/Test4/config/auspice_config.json"

try:
    with open(json_file) as ifile:
        config = json.load(ifile)
except json.decoder.JSONDecodeError as err:
    print("\tError message: '{}'".format(err.msg))
    print("\tLine number: '{}'".format(err.lineno))
    print("\tColumn number: '{}'".format(err.colno))

