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
	pairs_file = args.pairs_file
	dataset = args.dataset
	run = args.run
	dag = args.dag
	cores = args.cores
	cluster = args.cluster
	until = args.until
	verbose = args.verbose
	list_params_changes = args.list_params_changes
	keep_incomplete = args.keep_incomplete
	summary = args.summary

	logger.info("Pipeline started using: " + pairs_file)
	
	# Define initial config parameters
	params = {
		"dataset": dataset,
		"pairs_file": pairs_file,
		"output": 'pipeline',
		"cluster": False,
	}
	
	# Read config file
	with open('kasel-config.yml', "r") as stream:
		try:
			config = yaml.safe_load(stream)
		except yaml.YAMLError as exc:
			logger.error("Unable to load config file: " + str(exc))

	# Fetch cluster command
	if cluster == True:
		try:
			params['cluster'] = config['cluster']
		except:
			logger.warning("No cluster config parameter found, defaulting to qsub")
			params['cluster'] = 'qsub'
		
		logger.info("Cluster command: " + params['cluster'])

	# Check if specific rules called
	if until == None:
		until = []
	else:
		until = [until]
	
	# Fetching current diractory
	thisdir = os.path.abspath(os.path.dirname(__file__))
	tempdir = ''

	snakefile = os.path.join(thisdir, 'workflow.pairs','Snakefile')
	
	if not os.path.exists(snakefile):
		logger.error("Cannot find snakefile at {}\n".format(snakefile))
		sys.exit(-1)

	status = snakemake.snakemake(snakefile, printshellcmds=verbose, printreason=verbose, quiet=False, forceall=False, force_incomplete=True,
									list_params_changes=list_params_changes, keep_incomplete=keep_incomplete, summary=summary, keepgoing=True, latency_wait=100,
									workdir=tempdir, config=params, cores=cores, nodes=cores, lock=False, dryrun=run, use_conda=True, 
									cluster=params['cluster'], conda_frontend="mamba", printdag=dag, until=until)
	
	if status == True:
		logger.info("Pipeline finished")
	else:
		logger.error("Pipeline finished with errors")

