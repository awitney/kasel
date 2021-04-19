import sys
import sqlite3
import pandas as pd
import re
import os
import snakemake
import yaml
from collections import defaultdict

import custom_logger as custom_logger

from loguru import logger

def run(parser, args):
	
	meta_file = args.meta_file
	DATASET = args.dataset
	run = args.run
	cores = args.cores

	logger.info("Pipeline started using: " + meta_file)

	thisdir = os.path.abspath(os.path.dirname(__file__))
	tempdir = ''

	params = {
		"dataset": DATASET,
		"meta_file": meta_file,
		"output": 'pipeline',
		"cluster": 'qsub -V -l h_rt=48:00:00 -l mem={resources.memory} -pe smp {threads}'
	}

	snakefile = os.path.join(thisdir, 'workflow','Snakefile')
	
	if not os.path.exists(snakefile):
		logger.error("Cannot find snakefile at {}\n".format(snakefile))
		sys.exit(-1)

	status = snakemake.snakemake(snakefile, printshellcmds=False, quiet=False, forceall=False, force_incomplete=True,
									workdir=tempdir, config=params, nodes=cores, lock=False, dryrun=run, use_conda=True, cluster=params['cluster'], conda_frontend="conda")

	logger.info("Pipeline finished: " + str(status))

