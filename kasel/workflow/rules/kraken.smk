
rule run_kraken:
	threads:
		5
	resources:
		memory = config['run_kraken']['memory']
	params:
		krakendb = KRAKENDB
	input:
		r1 = lambda wildcards: get_seq(wildcards, 'forward'),
		r2 = lambda wildcards: get_seq(wildcards, 'reverse'),
	output:
		report = join(KRAKEN, '{sample}.report.txt'),
#		output = join(KRAKEN, '{sample}.output.txt')
	log:
		join(LOGS, 'kraken.{sample}.log')
	conda:
		"../envs/kraken.yml"
	shell:
		"""
		kraken2 --db {params.krakendb}  --paired --gzip-compressed --thread {threads} --use-names --report {output.report} --output - {input.r1} {input.r2} &> {log}
		"""
#		kraken2 --db {params.krakendb}  --paired --gzip-compressed --thread {threads} --use-names --confidence 0.1 --report {output.report} --output {output.output} {input.r1} {input.r2} &> {log}


rule run_kraken_extract:
	threads:
		5
	params:
		taxid = '{taxid}',
		r1 = join(DATA, DATASET + '-{taxid}', '{sample}_1.fastq'),
		r2 = join(DATA, DATASET + '-{taxid}', '{sample}_2.fastq'),
	input:
		r1 = lambda wildcards: get_seq(wildcards, 'forward'),
		r2 = lambda wildcards: get_seq(wildcards, 'reverse'),
		report = join(KRAKEN, '{sample}.report.txt'),
		output = join(KRAKEN, '{sample}.output.txt'),
	output:
		r1 = join(DATA, DATASET + '-{taxid}', '{sample}_1.fastq.gz'),
		r2 = join(DATA, DATASET + '-{taxid}', '{sample}_2.fastq.gz'),
	log:
		join(LOGS, 'kraken-extract-{taxid}.{sample}.log')
	conda:
		"../envs/kraken.yml"
	shell:
		"""
		extract_kraken_reads.py -k {input.output} --report {input.report} -s1 {input.r1} -s2 {input.r2} -o {params.r1} -o2 {params.r2} -t {params.taxid} --include-children --fastq-output
		gzip -S .gz {params.r1}
		gzip -S .gz {params.r2}
		"""
