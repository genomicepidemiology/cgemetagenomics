import os
import sys
import subprocess
import csv

from cgemetagenomics import kma

def metagenomics_pipeline(args):
    print("Starting the metagenomics pipeline...")

    # Check if output folder already exists
    output_dir = '/var/lib/cge/results/{}'.format(args.name)
    if os.path.exists(output_dir):
        sys.exit(
            f"Error: Output directory '{output_dir}' already exists. Please choose a different name or delete the existing directory.")

    if args.db_dir is None:
        if not os.path.exists('/var/lib/cge/database/cge_db'):
            sys.exit('Please install the cge_db. It should be located in /var/lib/cge/database/cge_db')
        else:
            args.db_dir = '/var/lib/cge/database/cge_db'
            print(f"Using CGE database directory: {args.db_dir}")

    if args.output is None:
        args.output = '/var/lib/cge/results/{}'.format(args.name)

    # Create output directory
    print(f"Creating output directory: {args.output}")
    os.system('mkdir -p ' + args.output)

    # Load pathogen species
    print("Loading pathogen species...")
    species = load_pathogen_species(args.db_dir + '/pathogen_strains.list')

    # Run KMA for bacteria alignment
    print("Running KMA for bacteria alignment...")
    kma.KMARunner(args.input,
                  args.output + "/bacteria_alignment",
                  args.db_dir + '/bac_db/bac_db',
                  "-ID 25 -md 1 -ont -1t1 -mem_mode -t 8").run()

    bacterial_results = read_tab_separated_file(args.output + "/bacteria_alignment.res")

    for hit in bacterial_results:
        if 'Escherichia coli' in hit['#Template']:
            e_coli_depth = find_max_depth_for_escherichia_coli(bacterial_results)
            print(f"Escherichia coli detected with depth {e_coli_depth}. Running virulence finder...")
            kma.KMARunner(args.input,
                          args.output + "/virulence",
                          args.db_dir + '/virulence_db/virulence_db',
                          "-ont -md {} -mem_mode -t 8".format(e_coli_depth / 2)).run()
            break

    # Run KMA for AMR
    print("Running KMA for AMR...")
    kma.KMARunner(args.input,
                  args.output + "/amr",
                  args.db_dir + '/resfinder_db/resfinder_db',
                  "-ont -md 3 -mem_mode -t 8").run()

    # Generate report
    print("Creating refined report...")
    report = create_refined_report(args.db_dir + '/phenotypes.txt', args.output, bacterial_results, species, args.name)
    with open(args.output + '/report.txt', 'w') as report_file:
        report_file.write(report)

    print("Metagenomics pipeline completed successfully. Report generated and stored in " + args.output + '/report.txt')
    return 'metagenomics_pipeline'


def merge_fastq_files_unix(source_directory, output_name):
    """
    Merge all fastq.gz files in the given directory using Unix commands and save the output with the specified name in the home directory.

    Args:
    source_directory (str): Path to the directory containing fastq.gz files.
    output_name (str): Name for the output file.
    """
    # Home directory path
    home_directory = os.path.expanduser('~')

    # Output file path with the specified name
    output_file = os.path.join(home_directory, f'{output_name}.fastq.gz')

    # Creating the Unix command for concatenation
    cmd = f'cat {source_directory}/*.fastq.gz > {output_file}'

    # Executing the command
    subprocess.run(cmd, shell=True, check=True)

    print(f"All files merged into {output_file}")

def create_refined_report(phenotype_file, output, bacterial_results, species, name):
    gene_data = read_tab_separated_file(phenotype_file)
    amr_results = read_tab_separated_file(output + '/amr.res')
    e_coli_found = any(extract_species(result.get('#Template', '')) == 'Escherichia coli' for result in bacterial_results)

    pathogens = []
    non_pathogens = []
    for result in bacterial_results:
        species_name = extract_species(result.get('#Template', ''))
        if species_name in species:
            pathogens.append(result)
        else:
            non_pathogens.append(result)

    phenotypes = set()
    for amr_result in amr_results:
        for gene in gene_data:
            if gene['Gene_accession no.'] == amr_result.get('#Template'):
                phenotypes.update(gene['Phenotype'].split(','))

    report = "Analysis Report: {}\n".format(name)
    report += "=" * 60 + "\n"

    # AMR Results Section
    report += "Antimicrobial Resistance (AMR) Findings:\n"
    report += "-" * 60 + "\n"
    for result in amr_results:
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
        for phenotype in sorted(phenotypes):
            report += f"• {phenotype.strip()}\n"
    else:
        report += "No phenotypes expected based on AMR genes.\n"
    report += "\n"

    # Virulence Factors for Escherichia coli Section
    if e_coli_found:
        virulence_results = read_tab_separated_file(output + '/virulence.res')
        report += "Virulence Factors for Escherichia coli:\n"
        report += "-" * 60 + "\n"
        for result in virulence_results:
            report += f"Template: {result['#Template']}\n"
            report += f"Identity: {result['Template_Identity'].strip()}, Coverage: {result['Template_Coverage'].strip()}, Depth: {result['Depth'].strip()}\n\n"
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