
rule run_bakta:
	threads:
		config['bakta']['threads']
	params:
		name = '{sample}',
		outdir = join(BAKTA, '{sample}'),
		db = BAKTADB
	input:
		fasta = join(SHOVILL, '{sample}', 'contigs.fa')
	output:
		fna = join(BAKTA, '{sample}', '{sample}.fna'),
		gff = join(BAKTA, '{sample}', '{sample}.gff3')
	log:
		join(LOGS, 'bakta.{sample}.log')
	conda:
		"../envs/bakta.yml"
	shell:
		"""
		bakta --threads {threads} --db {params.db} --prefix {params.name} --output {params.outdir} {input.fasta} --force &> {log}
		"""
		
#		bakta --threads {threads} --db {params.db} --locus-tag {params.name} --prefix {params.name} --output {params.outdir} {input.fasta} --force &> {log}
