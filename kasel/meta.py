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

	with open("config.yml", "r") as ymlfile:
		config = yaml.load(ymlfile)

	config["dataset"] = DATASET
	config["meta_file"] = meta_file
	config["output"] = 'pipeline'

	print("META", config)
	snakefile = os.path.join(thisdir, 'workflow','Snakefile')
	
	if not os.path.exists(snakefile):
		logger.error("Cannot find snakefile at {}\n".format(snakefile))
		sys.exit(-1)

#	logger2 = custom_logger.Logger()
#	logger2.quiet = True
	
	cluster_config = 'config.yml'

	cluster = 'qsub -V -l h_rt=48:00:00 -l mem={resources.memory} -pe smp {threads} -wd /scratch/scratch/rekgawi/'
	status = snakemake.snakemake(snakefile, printshellcmds=False, quiet=False, forceall=False, force_incomplete=True,
									workdir=tempdir, nodes=cores, lock=False, dryrun=run, use_conda=True, cluster=cluster, conda_frontend="conda")
#									workdir=tempdir, config=config, cores=cores, lock=False, log_handler=logger2.log_handler, dryrun=run, use_conda=True)


	logger.info("Pipeline finished: " + str(status))

