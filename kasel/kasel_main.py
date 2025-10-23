#!/usr/bin/env python

# Written by Adam Witney
# Thanks to Aaron Quinlan for the argparse implementation from poretools.

import argparse
import sys

from . import version

from loguru import logger

def run_subtool(parser, args):
	if args.command == 'meta':
		from . import meta as submodule
	elif args.command == 'sample':
		from . import sample as submodule
	elif args.command == 'pairs':
		from . import pairs as submodule
	# run the chosen submodule.
	submodule.run(parser, args)


class ArgumentParserWithDefaults(argparse.ArgumentParser):
	def __init__(self, *args, **kwargs):
		super(ArgumentParserWithDefaults, self).__init__(*args, **kwargs)
		self.add_argument("-q", "--quiet", help="Do not output warnings to stderr", action="store_true", dest="quiet")
		self.add_argument('-v', '--verbose', action="store_true", dest='verbose', help='output extra logging information')
		self.add_argument('--debug', action="store_true", dest='debug', help='output debugging information')
		self.add_argument('-t', '--cpus', type=int, dest='cpus', default=1, help='number of cpus')
		self.add_argument('-r', '--run', action="store_false", required=False, dest='run', help='Run pipeline')
		self.add_argument('--summary', action="store_true", required=False, dest='summary', help='Show summary output')

def init_pipeline_parser():
	"""Wraps the argparse parser initialisation.

	Returns
	-------
	argparse.ArgumentParser
		The initialised argparse Argument Parser for the pipeline
	"""

	parser = argparse.ArgumentParser(
	prog='kasel', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	subparsers = parser.add_subparsers(
		title='[sub-commands]', dest='command', parser_class=ArgumentParserWithDefaults)
	parser.add_argument("--version", help="Installed Artic version", action="version", version="%(prog)s " + str(version.__version__))

	# meta
	parser_meta = subparsers.add_parser('meta', help='Meta data file')

	parser_meta.add_argument(
		'-m', '--meta', required=True, dest='meta_file', help='Meta data file')
	parser_meta.add_argument(
		'-d', '--dataset', required=True, dest='dataset', help='dataset')
	parser_meta.add_argument(
		'--dag', action="store_true", required=False, dest='dag', help='Generate DAG diagram')
	parser_meta.add_argument(
		'-cl', '--cluster', action="store_true", required=False, dest='cluster', help='Run on cluster with qsub')
	parser_meta.add_argument(
		'-c', '--cores', required=False, type=int, dest='cores', default=1, help='Number of cores')
	parser_meta.add_argument(
		'-u', '--until', required=False, dest='until', help='define specific rules to run')
	parser_meta.set_defaults(func=run_subtool)

	# sample
	parser_sample = subparsers.add_parser('sample', help='Sample subparser')
	parser_sample.add_argument("--version", help="Installed Artic version", action="version", version="%(prog)s " + str(version.__version__))

	parser_sample.add_argument(
		'-s', '--samples', required=True, dest='samples_file', help='Sample list file')
	parser_sample.add_argument(
		'-d', '--dataset', required=True, dest='dataset', help='dataset')
	parser_sample.set_defaults(func=run_subtool)
	parser_sample.add_argument(
		'--dag', action="store_true", required=False, dest='dag', help='Generate DAG diagram')
	parser_sample.add_argument(
		'-cl', '--cluster', action="store_true", required=False, dest='cluster', help='Run on cluster with qsub')
	parser_sample.add_argument(
		'-c', '--cores', required=False, type=int, dest='cores', default=1, help='Number of cores')
	parser_sample.add_argument(
		'-u', '--until', required=False, dest='until', help='define specific rules to run')
	parser_sample.add_argument(
		'--nolegacy', required=False, action="store_true", dest='nolegacy', help='Switch off legacy pipeline')
	parser_sample.add_argument(
		'--nocallers', required=False, action="store_true", dest='nocallers', help='Switch off external callers')
	parser_sample.add_argument(
		'--list-params-changes', required=False, action="store_true", dest='list_params_changes', help='List params changes')
	parser_sample.add_argument(
		'--keep-incomplete', required=False, action="store_true", dest='keep_incomplete', help='Keep incomplete files')
	parser_sample.add_argument(
		'--fastlin', required=False, dest='fastlin', help='Data directory for fastlin tool')
	parser_sample.add_argument(
		'--do_assembly', required=False, action="store_true", dest='do_assembly', help='Perform Shovill assembly and Bakta annotation')
	parser_sample.add_argument(
		'--conda_frontend', required=False, type=str, dest='conda_frontend', default="conda", help='Specify if using conda or mamba')

	parser_sample.set_defaults(func=run_subtool)

	# pairs
	parser_pairs = subparsers.add_parser('pairs', help='Sample subparser')
	parser_pairs.add_argument("--version", help="Installed Artic version", action="version", version="%(prog)s " + str(version.__version_pairs__))

	parser_pairs.add_argument(
		'-s', '--samples', required=True, dest='samples_file', help='Sample list file')
	parser_pairs.add_argument(
		'-p', '--pairs', required=True, dest='pairs_file', help='Sample list file')
	parser_pairs.add_argument(
		'-d', '--dataset', required=True, dest='dataset', help='dataset')
	parser_pairs.set_defaults(func=run_subtool)
	parser_pairs.add_argument(
		'--dag', action="store_true", required=False, dest='dag', help='Generate DAG diagram')
	parser_pairs.add_argument(
		'-cl', '--cluster', action="store_true", required=False, dest='cluster', help='Run on cluster with qsub')
	parser_pairs.add_argument(
		'-c', '--cores', required=False, type=int, dest='cores', default=1, help='Number of cores')
	parser_pairs.add_argument(
		'-u', '--until', required=False, dest='until', help='define specific rules to run')
	parser_pairs.add_argument(
		'--list-params-changes', required=False, action="store_true", dest='list_params_changes', help='List params changes')
	parser_pairs.add_argument(
		'--keep-incomplete', required=False, action="store_true", dest='keep_incomplete', help='Keep incomplete files')

	# return the parser
	return parser


def main():

	# init the pipeline parser
	parser = init_pipeline_parser()

	# collect the args
	args = parser.parse_args(sys.argv[1:])

	# set the logging level
	log_level = 'INFO'

	if args.debug:
		log_level = 'DEBUG'
	
	if args.quiet:	
		log_level = 'ERROR'

	logger.remove()
	logger.add(sys.stderr, level=log_level)

# if args.quiet:
#	logger.setLevel(logging.ERROR)

	# run the subcommand or print usage if no subcommand provided
	if args.command:
		args.func(parser, args)
	else:
		parser.print_usage()


if __name__ == "__main__":
	main()
