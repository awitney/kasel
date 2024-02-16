
rule tbprofiler_pe:
	threads:
		config['tbprofiler']['threads']
	resources:
		memory = config['tbprofiler']['memory']
	params:
		prefix = '{sample}',
		outdir = join(CALLERS, 'TB-profiler', '{sample}')
	input:
#		r1 = join(DATA, '{sample}_1.fastq.gz'),
#		r2 = join(DATA, '{sample}_2.fastq.gz')
		r1 = lambda wildcards: get_seq(wildcards, 'forward'),
		r2 = lambda wildcards: get_seq(wildcards, 'reverse'),
	output:
		txt1 = join(CALLERS, 'TB-profiler', '{sample}', 'results', '{sample}.results.json'),
		txt2 = join(CALLERS, 'TB-profiler', 'results', '{sample}.results.json'),
		txt = join(CALLERS, 'TB-profiler', '{sample}', 'results', '{sample}.results.txt'),
	log:
		join(LOGS, 'tbprofiler_pe.{sample}.log')
	message:
		"Running tbprofiler_pe for {wildcards.sample}"
	conda:
		"../envs/tb-profiler.yml"
	shell:
		"""
		tb-profiler profile --threads {threads} -1 {input.r1} -2 {input.r2}  -p {params.prefix} --txt --dir {params.outdir} &> {log}
		cp {output.txt1} {output.txt2}
		"""

rule tbprofiler_collate:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	params:
		indir = join(CALLERS, 'TB-profiler', 'results'),
		prefix = join(CALLERS, 'TB-profiler', 'tbprofiler'),
	input:
		expand(join(CALLERS, 'TB-profiler', 'results', '{sample}.results.json'), sample=SAMPLES),
	output:
		join(CALLERS, 'TB-profiler', 'tbprofiler.txt')
	log:
		join(LOGS, 'tbprofiler_collate.log')
	message:
		"Running tbprofiler_collate"
	conda:
		"../envs/tb-profiler.yml"
	shell:
		"""
		tb-profiler collate --dir {params.indir} --prefix {params.prefix} &> {log}
		"""

