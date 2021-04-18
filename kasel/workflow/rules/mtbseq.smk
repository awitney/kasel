
rule mtbseq:
	threads:
		config['mtbseq']['threads']
	resources:
		memory = config['mtbseq']['memory']
	params:
		data   = join(DATA, DATASET),
		outdir = join(CALLERS, DATASET, 'MTBSeq'),
	input:
		reads = READS1
	output:
		join(CALLERS, DATASET, 'MTBSeq', 'Classification/Strain_Classification.tab'),
	conda:
		"../envs/mtbseq.yml"
	log:
		"logs/mtbseq." + DATASET + ".log"
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
		data   = join(DATA, DATASET),
		outdir = join(CALLERS, DATASET, 'MTBSeq'),
	input:
	output:
		join(CALLERS, DATASET, 'MTBSeq', 'samples.txt'),
	run:
		"""
		print("We are printing samples.txt here", SAMPLES)
		"""
