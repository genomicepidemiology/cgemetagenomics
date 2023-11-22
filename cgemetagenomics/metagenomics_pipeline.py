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

    for hit in bacterial_results:
        if 'Escherichia coli' in hit['#Template']:
            e_coli_depth = find_max_depth_for_escherichia_coli(bacterial_results)
            # run virulence finder
            kma.KMARunner(args.input,
                          args.output + "/virulence",
                          args.db_dir + '/virulence_db/virulence_db',
                          "-ont -md {} -mem_mode -t 8".format(e_coli_depth / 2)).run()
            break

    kma.KMARunner(args.input,
                  args.output + "/amr",
                  args.db_dir + '/resfinder_db/resfinder_db',
                  "-ont -md 3 -mem_mode -t 8").run()

    amr_results = read_tab_separated_file(args.output + "/amr.res")

    report = create_refined_report(args.db_dir + '/phenotypes.txt', args.output, bacterial_results)
    print (report)

    #Parse bacterial alignment and output those above a set of thresholds

    return 'isolate_pipeline'


def create_refined_report(amr_file_path, bacterial_results, species, virulence_file_path=None):
    gene_data = read_tab_separated_file(amr_file_path)

    # Determine pathogens based on species list
    pathogens = []
    non_pathogens = []
    for result in bacterial_results:
        if extract_species(result['#Template']) in species:
            pathogens.append(result)
        else:
            non_pathogens.append(result)

    phenotypes = set()
    for amr_result in bacterial_results:
        for gene in gene_data:
            if gene['Gene_accession no.'] == amr_result['#Template']:
                phenotypes.update(gene['Phenotype'].split(','))

    report = "Sample Analysis Report\n"
    report += "=" * 60 + "\n"

    # AMR Results Section
    report += "Antimicrobial Resistance (AMR) Findings:\n"
    report += "-" * 60 + "\n"
    for result in bacterial_results:
        report += f"Template: {result['#Template']}\n"
        report += f"Identity: {result['Template_Identity'].strip()}, Coverage: {result['Template_Coverage'].strip()}, Depth: {result['Depth'].strip()}\n\n"

    # Pathogens Found Section
    report += "Identified Pathogens:\n"
    report += "-" * 60 + "\n"
    for pathogen in pathogens:
        report += f"Template: {pathogen['#Template']}\n"
        report += f"Depth: {pathogen['Depth'].strip()}, Coverage: {pathogen['Template_Coverage'].strip()}, Identity: {pathogen['Template_Identity'].strip()}, Length: {pathogen['Template_length'].strip()}\n\n"

    # Non-Pathogens Section
    report += "Non-Pathogenic Bacteria Detected:\n"
    report += "-" * 60 + "\n"
    for non_pathogen in non_pathogens:
        report += f"Template: {non_pathogen['#Template']}\n"
        report += f"Depth: {non_pathogen['Depth'].strip()}, Coverage: {non_pathogen['Template_Coverage'].strip()}, Identity: {non_pathogen['Template_Identity'].strip()}, Length: {non_pathogen['Template_length'].strip()}\n\n"

    # Expected Phenotypes Based on AMR Genes Section
    report += "Expected Phenotypes Based on AMR Genes:\n"
    report += "-" * 60 + "\n"
    if phenotypes:
        report += ', '.join(phenotypes) + "\n"
    else:
        report += "No phenotypes expected based on AMR genes.\n"
    report += "\n"

    # Virulence Factors for Escherichia coli Section
    if 'Escherichia coli' in species and virulence_file_path:
        virulence_results = read_tab_separated_file(virulence_file_path)
        report += "Virulence Factors for Escherichia coli:\n"
        report += "-" * 60 + "\n"
        for result in virulence_results:
            report += str(result) + "\n"
    else:
        report += "No virulence factors analysis for Escherichia coli.\n"

    return report

def find_max_depth_for_escherichia_coli(bacterial_results):
    max_depth = 0.0

    for hit in bacterial_results:
        if 'Escherichia coli' in hit['#Template']:
            # Convert depth to a float and compare with the current max
            depth = float(hit['Depth'].strip())
            if depth > max_depth:
                max_depth = depth

    return max_depth

def load_tsv(file_path):
    with open(file_path, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        return list(reader)


def load_pathogen_species(strain_file):
    strains = []
    with open(strain_file, 'r') as file:
        for line in file:
            line = line.strip().split('\t')
            species = line[1].split(' ')[0] + ' ' + line[1].split(' ')[1]
            strains.append(species)
    return strains


def read_tab_separated_file(file_path):
    with open(file_path, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        return list(reader)

def extract_species(template_str):
    return ' '.join(template_str.split()[1:3])