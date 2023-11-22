import os
import sys
import subprocess
import csv

from cgemetagenomics import kma

def metagenomics_pipeline(args):
    os.system('mkdir ' + args.output)
    # Check if kma is installed
    species = load_pathogen_species(args.db_dir + '/pathogen_strains.list')
    kma.KMARunner(args.input,
              args.output + "/bacteria_alignment",
              args.db_dir + '/bac_db/bac_db',
              "-ID 25 -md 1 -ont -1t1 -mem_mode -t 8").run()

    baterial_results = read_tab_separated_file(args.output + "/bacteria_alignment.res")

    filter_and_print_hits(baterial_results, species)

    sys.exit()

    kma.KMARunner(args.input,
                  args.output + "/amr",
                  args.db_dir + '/resfinder_db/resfinder_db',
                  "-ont -md 3 -mem_mode -t 8").run()

    amr_results = read_tab_separated_file(args.output + "/amr.res")


    #Parse bacterial alignment and output those above a set of thresholds

    return 'isolate_pipeline'

def filter_and_print_hits(bacterial_results, species_list):
    for hit in bacterial_results:
        # Extract species name from the '#Template' field
        template = hit['#Template']
        species_name = ' '.join(template.split()[1:3])
        print (species_name)

        # Check if the species name is in the list of species
        if species_name in species_list:
            print (hit, 'found')

def load_pathogen_species(strain_file):
    strains = []
    with open(strain_file, 'r') as file:
        for line in file:
            line = line.strip().split('\t')
            species = line[1].split(' ')[0] + ' ' + line[1].split(' ')[1]
            strains.append(species)
    return strains


def read_tab_separated_file(file_path):
    results = []
    with open(file_path, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            results.append(row)
    return results