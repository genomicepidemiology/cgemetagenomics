# cgemetagenomics

# Installation

The following Conda channels are required:

`conda config --add channels defaults`

`conda config --add channels bioconda`

`conda config --add channels conda-forge`

Mamba can be installed with:

`conda install -c conda-forge mamba`

For a fast install of cgemetagenomics, use mamba:
`mamba install -c genomicepidemiology cgemetagenomics`


# Datebase

Download the cge_db database:

`wget https://cge.food.dtu.dk/services/MINTyper/cge_db.tar.gz`

`tar -xvzf cge_db.tar.gz`

`sudo mkdir -m 777 /var/lib/cge`

`sudo mkdir -m 777 /var/lib/cge/database`

`mv cge_db /var/lib/cge/database/cge_db`



# Usage

Standard Usage:
`cgemetagenomics -i <input_file> -o <output_file>`

If your cge_db is not stored in /opt/cge/cge_db:
`cgemetagenomics -i <input_file> -o <output_file> -db_dir <path_to_cge_db>`

Help Message:
`cgemetagenomics -h`
