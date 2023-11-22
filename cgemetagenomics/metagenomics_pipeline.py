import os
import sys
import subprocess
import csv

from cgemetagenomics import kma

def metagenomics_pipeline(args):
    os.system('mkdir ' + args.output)
    # Check if kma is installed

    load_pathogen_strains(args.db_dir + '/pathogen_strains.list')
    sys.exit()
    kma.KMARunner(args.input,
              args.output + "/bacteria_alignment",
              args.db_dir + '/bac_db/bac_db',
              "-ID 25 -md 1 -ont -1t1 -mem_mode -t 8").run()

    baterial_results = read_tab_separated_file(args.output + "/bacteria_alignment.res")

    kma.KMARunner(args.input,
                  args.output + "/amr",
                  args.db_dir + '/resfinder_db/resfinder_db',
                  "-ont -md 3 -mem_mode -t 8").run()

    amr_results = read_tab_separated_file(args.output + "/amr.res")

    #Parse bacterial alignment and output those above a set of thresholds

    return 'isolate_pipeline'

def load_pathogen_strains(strain_file):
    strains = []
    with open(strain_file, 'r') as file:
        for line in file:
            line = line.strip().split('\t')
            print (line)

def read_tab_separated_file(file_path):
    results = []
    with open(file_path, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            results.append(row)
    return results