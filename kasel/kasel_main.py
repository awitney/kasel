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

	# run the chosen submodule.
	submodule.run(parser, args)


class ArgumentParserWithDefaults(argparse.ArgumentParser):
	def __init__(self, *args, **kwargs):
		super(ArgumentParserWithDefaults, self).__init__(*args, **kwargs)
		self.add_argument("-q", "--quiet", help="Do not output warnings to stderr",
						action="store_true",
						dest="quiet")
		self.add_argument('-v', '--verbose', action="store_true", dest='verbose',
						help='output extra logging information')
		self.add_argument('--debug', action="store_true", dest='debug',
						help='output debugging information')
		self.add_argument('-t', '--cpus', type=int, dest='cpus', default=1,
						help='number of cpus')


def init_pipeline_parser():
	"""Wraps the argparse parser initialisation.

	Returns
	-------
	argparse.ArgumentParser
		The initialised argparse Argument Parser for the pipeline
	"""

	parser = argparse.ArgumentParser(
	prog='kasel', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-v", "--version", help="Installed Artic version",
		action="version",
		version="%(prog)s " + str(version.__version__))
	subparsers = parser.add_subparsers(
		title='[sub-commands]', dest='command', parser_class=ArgumentParserWithDefaults)

	# meta
	parser_meta = subparsers.add_parser('meta', help='Meta data file')

	parser_meta.add_argument(
		'-m', '--meta', required=True, dest='meta_file', help='Meta data file')
	parser_meta.add_argument(
		'-d', '--dataset', required=True, dest='dataset', help='dataset')
	parser_meta.add_argument(
		'-r', '--run', action="store_false", required=False, dest='run', help='Run pipeline')
	parser_meta.add_argument(
		'--dag', action="store_true", required=False, dest='dag', help='Generate DAG diagram')
	parser_meta.add_argument(
		'-cl', '--cluster', action="store_true", required=False, dest='cluster', help='Run on cluster with qsub')
	parser_meta.add_argument(
		'-c', '--cores', required=False, type=int, dest='cores', default=1, help='Number of cores')
	parser_meta.add_argument(
		'-u', '--until', required=False, dest='until', help='define specific rules to run')
	parser_meta.set_defaults(func=run_subtool)


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
