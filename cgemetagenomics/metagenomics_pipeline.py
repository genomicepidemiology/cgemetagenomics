import os
import sys
import subprocess

from cgemetagenomics import kma

def isolate_pipeline(args):
    os.system('mkdir ' + args.output)
    # Check if kma is installed
    kma.KMARunner(args.input,
              args.output + "/bacteria_alignment",
              args.db_dir + '/bac_db/bac_db',
              "-ID 25 -ont -1t1 -mem_mode -t 8").run()

    #Parse bacterial alignment and output those above a set of thresholds

    return 'isolate_pipeline'