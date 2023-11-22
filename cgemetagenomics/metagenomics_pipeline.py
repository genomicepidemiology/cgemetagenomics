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

    bacterial_results = read_tab_separated_file(args.output + "/bacteria_alignment.res")

    pathogens_found = filter_and_print_hits(bacterial_results, species)

    if 'Escherichia coli' in pathogens_found:

        e_coli_depth = find_max_depth_for_escherichia_coli(bacterial_results)
        #run virulence finder
        kma.KMARunner(args.input,
                      args.output + "/virulence",
                      args.db_dir + '/virulence_db/virulence_db',
                      "-ont -md {} -mem_mode -t 8".format(e_coli_depth/2)).run()

    kma.KMARunner(args.input,
                  args.output + "/amr",
                  args.db_dir + '/resfinder_db/resfinder_db',
                  "-ont -md 3 -mem_mode -t 8").run()

    amr_results = read_tab_separated_file(args.output + "/amr.res")

    report = create_refined_report(amr_results, bacterial_results, pathogens_found)
    print (report)

    #Parse bacterial alignment and output those above a set of thresholds

    return 'isolate_pipeline'

def create_refined_report(amr_results, pathogens_found, non_pathogens):
    report = "Sample Analysis Report\n"
    report += "=" * 60 + "\n"

    # AMR Results Section
    report += "Antimicrobial Resistance (AMR) Findings:\n"
    report += "-" * 60 + "\n"
    if amr_results:
        for result in amr_results:
            report += f"Template: {result['#Template']}\n"
            report += f"Identity: {result['Template_Identity'].strip()}, Coverage: {result['Template_Coverage'].strip()}, Depth: {result['Depth'].strip()}\n\n"
    else:
        report += "No AMR genes detected.\n"
    report += "\n"

    # Pathogens Found Section
    report += "Identified Pathogens:\n"
    report += "-" * 60 + "\n"
    if pathogens_found:
        for pathogen in pathogens_found:
            report += f"Template: {pathogen['#Template']}\n"
            report += f"Depth: {pathogen['Depth'].strip()}, Coverage: {pathogen['Template_Coverage'].strip()}, Identity: {pathogen['Template_Identity'].strip()}, Length: {pathogen['Template_length'].strip()}\n\n"
    else:
        report += "No pathogens detected.\n"
    report += "\n"

    # Non-Pathogens Section
    report += "Non-Pathogenic Bacteria Detected:\n"
    report += "-" * 60 + "\n"
    if non_pathogens:
        for non_pathogen in non_pathogens:
            report += f"Template: {non_pathogen['#Template']}\n"
            report += f"Depth: {non_pathogen['Depth'].strip()}, Coverage: {non_pathogen['Template_Coverage'].strip()}, Identity: {non_pathogen['Template_Identity'].strip()}, Length: {non_pathogen['Template_length'].strip()}\n\n"
    else:
        report += "No non-pathogenic bacteria detected.\n"

    return reportdef find_max_depth_for_escherichia_coli(bacterial_results):
    max_depth = 0.0

    for hit in bacterial_results:
        if 'Escherichia coli' in hit['#Template']:
            # Convert depth to a float and compare with the current max
            depth = float(hit['Depth'].strip())
            if depth > max_depth:
                max_depth = depth

    return max_depth

def filter_and_print_hits(bacterial_results, species_list):
    pathogens_found = []
    for hit in bacterial_results:
        # Extract species name from the '#Template' field
        template = hit['#Template']
        species_name = ' '.join(template.split()[1:3])

        # Check if the species name is in the list of species
        if species_name in species_list:
            pathogens_found.append(species_name)
    return pathogens_found

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