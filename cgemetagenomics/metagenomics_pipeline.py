import os
import sys
import subprocess
import csv

from cgemetagenomics import kma

def metagenomics_pipeline(args):
    os.system('mkdir ' + args.output)
    # Check if kma is installed
    kma.KMARunner(args.input,
              args.output + "/bacteria_alignment",
              args.db_dir + '/bac_db/bac_db',
              "-ID 25 -md 1 -ont -1t1 -mem_mode -t 8").run()

    baterial_results = read_tab_separated_file(args.output + "/bacteria_alignment.res")
    for item in baterial_results:
        print (item)
    #Parse bacterial alignment and output those above a set of thresholds

    return 'isolate_pipeline'

def read_tab_separated_file(file_path):
    results = []
    with open(file_path, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            results.append(row)
    return results