
rule tbprofiler_pe:
	threads:
		config['tbprofiler']['threads']
	resources:
		memory = config['tbprofiler']['memory']
	params:
		prefix = '{sample}',
		outdir = join(CALLERS, DATASET, 'TB-profiler')
	input:
#		r1 = join(DATA, DATASET, '{sample}_1.fastq.gz'),
#		r2 = join(DATA, DATASET, '{sample}_2.fastq.gz')
		r1 = lambda wildcards: get_seq(wildcards, 'forward'),
		r2 = lambda wildcards: get_seq(wildcards, 'reverse'),
	output:
		join(CALLERS, DATASET, 'TB-profiler', 'results', '{sample}.results.txt')
	log:
		join(LOGS, DATASET, 'tbprofiler.{sample}.log')
	conda:
		"../envs/tb-profiler.yml"
	shell:
		"""
		tb-profiler profile --threads {threads} -1 {input.r1} -2 {input.r2}  -p {params.prefix} --txt --dir {params.outdir} 2>&1 > {log}
		"""

rule tbprofiler_collate:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	params:
		indir = join(CALLERS, DATASET, 'TB-profiler', 'results'),
		prefix = join(CALLERS, DATASET, 'TB-profiler', 'tbprofiler'),
	input:
		expand(join(CALLERS, DATASET, 'TB-profiler', 'results', '{sample}.results.txt'), sample=SAMPLES),
	output:
		join(CALLERS, DATASET, 'TB-profiler', 'tbprofiler.txt')
	conda:
		"../envs/tb-profiler.yml"
	shell:
		"""
		tb-profiler collate --dir {params.indir} --prefix {params.prefix}
		"""

