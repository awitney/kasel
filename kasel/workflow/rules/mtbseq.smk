
rule mtbseq:
	threads:
		config['mtbseq']['threads']
	resources:
		memory = config['mtbseq']['memory']
	params:
		data   = DATA,
		outdir = join(CALLERS, 'MTBSeq'),
	input:
		reads = READS1
	output:
		join(CALLERS, 'MTBSeq', 'Classification/Strain_Classification.tab'),
	conda:
		"../envs/mtbseq.yml"
	log:
		join(LOGS, "mtbseq.log")
	shell:
		"""
		cd {params.outdir}
		ln -s ../../../../{params.data}/*.gz ./
		rename _1.fastq.gz _lib1_R1.fastq.gz *_1.fastq.gz
		rename _2.fastq.gz _lib1_R2.fastq.gz *_2.fastq.gz
		MTBseq --threads {threads} --project MTBseq --step TBfull --continue >> ../../../../{log} 2>&1
		"""

rule mtbseq_samples:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	params:
		data   = DATA,
		outdir = join(CALLERS, 'MTBSeq'),
	input:
	output:
		join(CALLERS, 'MTBSeq', 'samples.txt'),
	run:
		"""
		print("We are printing samples.txt here", SAMPLES)
		"""
