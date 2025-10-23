
rule run_shovill:
	threads:
		config['shovill']['threads']
	params:
		outdir = join(SHOVILL, '{sample}')
	input:
		r1 = lambda wildcards: get_seq(wildcards, 'forward'),
		r2 = lambda wildcards: get_seq(wildcards, 'reverse'),
	output:
		fasta = join(SHOVILL, '{sample}/contigs.fa')
	conda:
		"../envs/shovill.yml"
	shell:
		"shovill --outdir {params.outdir} --R1 {input.r1} --R2 {input.r2} --cpus {threads} --force --ram 32"

