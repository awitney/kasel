
rule run_fastin:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	params:
		datadir = config['fastlin'],
		barcodes = config['fastlin_barcodes']
	output:
		out = join(RESULTS, 'fastlin.txt'),
		versions = join(RESULTS, 'fastlin.versions.txt')
	log:
		join(LOGS, 'fastlin.log')
	message:
		"Running fastlin"
	conda:
		"../envs/fastlin.yml"
	shell:
		"""
		fastlin -d {params.datadir} -b {params.barcodes} -o {output.out} &> {log}
		fastlin --version > {output.versions}
		"""
