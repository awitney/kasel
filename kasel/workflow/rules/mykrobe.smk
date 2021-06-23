
rule mykrobe_pe:
	threads:
		1
	params:
		sample = '{sample}',
		skel_dir = join(CALLERS, DATASET, 'mykrobe', '{sample}')
	input:
#		r1 = join(DATA, DATASET, '{sample}_1.fastq.gz'),
#		r2 = join(DATA, DATASET, '{sample}_2.fastq.gz')
		r1 = lambda wildcards: get_seq(wildcards, 'forward'),
		r2 = lambda wildcards: get_seq(wildcards, 'reverse'),
	output:
		join(CALLERS, DATASET, 'mykrobe', '{sample}.results.csv')
	conda:
		"../envs/mykrobe.yml"
	shell:
		"""
		mykrobe predict --skeleton_dir {params.skel_dir} {params.sample} tb -1 {input.r1} {input.r2} > {output}
		"""

rule mykrobe_se:
	threads:
		1
	params:
		sample = '{sample}',
		skel_dir = join(CALLERS, DATASET, 'mykrobe', '{sample}')
	input:
		r1 = join(DATA, DATASET, '{sample}.fastq.gz'),
	output:
#		join(CALLERS, DATASET, 'mykrobe', '{sample}.results.csv')
	conda:
		"../envs/mykrobe.yml"
	shell:
		"""
		mykrobe predict --skeleton_dir {params.skel_dir} {params.sample} tb -1 {input.r1} > {output}
		"""
