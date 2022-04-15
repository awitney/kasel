import sys
import pandas as pd
import re
import os
import snakemake
import yaml
from collections import defaultdict

#import custom_logger as custom_logger

from loguru import logger

def run(parser, args):
	
	samples_file = args.samples_file
	DATASET = args.dataset
	run = args.run
	dag = args.dag
	cores = args.cores
	cluster = args.cluster
	until = args.until
	nolegacy = args.nolegacy
	verbose = args.verbose
	list_params_changes = args.list_params_changes
	summary = args.summary

	logger.info("Pipeline started using: " + samples_file)

	thisdir = os.path.abspath(os.path.dirname(__file__))
	tempdir = ''

	params = {
		"dataset": DATASET,
		"samples_file": samples_file,
		"output": 'pipeline',
		"cluster": False,
		"nolegacy": nolegacy
	}
	
	if cluster == True:
		params['cluster'] = "qsub -V -l h_rt=48:00:00 -l mem={resources.memory} -pe smp {threads}"
	
	if until == None:
		until = []
	else:
		until = [until]
	
	snakefile = os.path.join(thisdir, 'workflow','Snakefile')
	
	if not os.path.exists(snakefile):
		logger.error("Cannot find snakefile at {}\n".format(snakefile))
		sys.exit(-1)

	status = snakemake.snakemake(snakefile, printshellcmds=verbose, printreason=verbose, quiet=False, forceall=False, force_incomplete=True,
									list_params_changes=list_params_changes, summary=summary, keepgoing=True,
									workdir=tempdir, config=params, cores=cores, nodes=cores, lock=False, dryrun=run, use_conda=True, 
									cluster=params['cluster'], conda_frontend="conda", printdag=dag, until=until)

	logger.info("Pipeline finished")

