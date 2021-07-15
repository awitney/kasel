
rule alignment_pe_old:
	threads:
		config['alignment_pe']['threads']
	resources:
		memory = config['alignment_pe']['memory']
	params:
		tempf = 'temp.{sample}'
	input:
		genome = join(GENOMES, '{ref}.fna'),
#		reads = lambda wildcards: READS[wildcards.sample],
#		r1 = join(DATA, DATASET, '{sample}_1.fastq.gz'),
#		r2 = join(DATA, DATASET, '{sample}_2.fastq.gz')
		r1 = lambda wildcards: get_seq(wildcards, 'forward'),
		r2 = lambda wildcards: get_seq(wildcards, 'reverse'),
	output:
		bam = join(ALIGNMENTS, DATASET, OLD, '{ref}_{sample}.bam')
	conda:
		"../envs/alignment-0.1.20.yml"
	shell:
		"bwa mem -t {threads} {input.genome} {input.r1} {input.r2} | samtools view -bS - | samtools sort -o - {params.tempf} - | samtools rmdup - - > {output.bam} && samtools index {output.bam}"

rule site_calling_old:
	threads:
		config['site_calling']['threads']
	resources:
		memory = config['site_calling']['memory']
	input:
		genome = join(GENOMES, '{ref}.fna'),
#		bam = rules.alignment.output.bam
		bam = join(ALIGNMENTS, DATASET, OLD, '{ref}_{sample}.bam'),
	output:
		vcf = join(VCF, DATASET, OLD, '{ref}_{sample}.all.vcf.gz')
	conda:
		"../envs/alignment-0.1.20.yml"
	shell:
		"""
		samtools mpileup -gf {input.genome} {input.bam} | bcftools view -cg - | bgzip > {output.vcf} && tabix -p vcf {output.vcf}
		"""

rule variant_calling_old:
	threads:
		config['variant_calling']['threads']
	resources:
		memory = config['variant_calling']['memory']
	params:
		ref = '{ref}.3',
		bcf = join(VCF, DATASET, OLD, '{ref}_{sample}.bcf')
	input:
		genome = join(GENOMES, '{ref}.fna'),
#		bam = rules.alignment.output.bam
		bam = join(ALIGNMENTS, DATASET, OLD, '{ref}_{sample}.bam'),
	output:
		vcf = join(VCF, DATASET, OLD, 'variants', '{ref}_{sample}.vcf.gz')
	conda:
		"../envs/alignment-0.1.20.yml"
	shell:
		"""
		samtools mpileup -ugf {input.genome} {input.bam} | bcftools view -bvcg - > {params.bcf} && bcftools view {params.bcf} | vcfutils.pl varFilter -D 20000 | sed 's/NC_000962.3/Chromosome/g' | bgzip > {output.vcf} && tabix -p vcf {output.vcf};
		rm {params.bcf}
		"""

rule snp_annotation_old:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	input:
		join(VCF, DATASET, OLD, 'variants', '{ref}_{sample}.vcf.gz'),
	output:
		join(VCF, DATASET, OLD, 'variants', 'annotated', '{ref}_{sample}.ann.vcf.gz'),
	conda:
		"../envs/alignment-0.1.20.yml"
	shell:
		"""
        snpEff eff -no-downstream -no-upstream -no-utr -o vcf Mycobacterium_tuberculosis_h37rv {input} | bgzip > {output} && bcftools index {output}
		"""

rule snp_report_all_old:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	input:
		tsv = expand(join(VCF, DATASET, OLD, 'variants/annotated', '{ref}_{sample}.tsv'), ref=REF, sample=SAMPLES),
	output:
		tsv = join(RESULTS, DATASET, OLD, 'variants', '{ref}_' + DATASET + '.tsv'),
	conda:
		"../envs/alignment.yml"
	shell:
		"""
		echo "Sample\tPos\tGT\tREF\tALT\tType\tQUAL\tFilter\tDP\tDP4\tANN\t" > {output.tsv}
		perl -p -e 's/.+\/NC_000962_(.+).bam/$1/' {input.tsv} >> {output.tsv}
		"""

rule snp_report_old:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	params:
		string = join(ALIGNMENTS, DATASET, OLD, '{ref}_')
	input:
		vcf = join(VCF, DATASET, OLD, 'variants/annotated', '{ref}_{sample}.ann.vcf.gz'),
	output:
		tsv = join(VCF, DATASET, OLD, 'variants/annotated', '{ref}_{sample}.tsv'),
	conda:
		"../envs/alignment.yml"
	shell:
		"""
		bcftools index -f {input.vcf}
		bcftools filter -sFAIL -g3 -G10 -e '%QUAL<30 || FORMAT/DP<4' {input.vcf} | \
			bcftools query -f '[%SAMPLE]\t%POS\t[%GT]\t%REF\t%ALT{{0}}\t%TYPE\t%QUAL\t%FILTER\t%INFO/DP\t[%INFO/DP4]\t[%INFO/ANN]\n' -i 'GT="alt"' - | \
			perl -p -e 's/\|/\t/g' > {output.tsv} 
		"""
#			perl -p -e 's/\|/\t/g' | perl -p -e 's|{params.string}(.+).bam|$1|' > {output.tsv} 

rule snp_report_resistance_all_old:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	input:
		tsv = expand(join(VCF, DATASET, OLD, 'variants/annotated/resistance', '{ref}_{{drug}}_{sample}.tsv'), ref=REF, sample=SAMPLES),
	output:
		tsv = join(RESULTS, DATASET, OLD, 'variants/resistance', '{ref}_' + DATASET + '_{drug}.tsv'),
	conda:
		"../envs/alignment.yml"
	shell:
		"""
		echo "Sample\tPos\tGT\tREF\tALT\tType\tQUAL\tFilter\tDP\tDP4\tANN\t" > {output.tsv}
		cat {input.tsv} >> {output.tsv}
		"""

rule snp_report_resistance_all_summary_old:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	params:
		finddir = join(VCF, DATASET, OLD, 'variants/annotated/resistance/'),
		string = join(VCF, DATASET, OLD, 'variants/annotated/resistance/' + REF + '_'),
	input:
		tsv1 = expand(join(RESULTS, DATASET, OLD, 'variants/resistance', '{ref}_' + DATASET + '_BDQ.tsv'), ref=REF),
		tsv2 = expand(join(RESULTS, DATASET, OLD, 'variants/resistance', '{ref}_' + DATASET + '_PTM.tsv'), ref=REF),
	output:
		out = join(RESULTS, DATASET, OLD, 'variants/resistance', '{ref}_' + DATASET + '.tsv'),
	conda:
		"../envs/alignment.yml"
	shell:
		"""
		wc -l  $(find {params.finddir} -name "*_BDQ*.tsv") | perl -p -e 's:[ ]+(\d+)[ ]+{params.string}(PTM|BDQ)_(.+).tsv:$3\t$2\t$1:' | sort -k2 > {output.out}
		wc -l  $(find {params.finddir} -name "*_PTM*.tsv") | perl -p -e 's:[ ]+(\d+)[ ]+{params.string}(PTM|BDQ)_(.+).tsv:$3\t$2\t$1:' | sort -k2 >> {output.out}
		"""

rule snp_report_resistance_old:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	params:
		string = join(ALIGNMENTS, DATASET, OLD, '{ref}_'),
		bed    = 'kasel/kasel/workflow/data/snps.Chromosome-{drug}.bed',
		tmp    = join(VCF, DATASET, OLD, 'variants', 'annotated', '{ref}_{sample}.ann.vcf.gz.tmp'),
	input:
		vcf = join(VCF, DATASET, OLD, 'variants', 'annotated', '{ref}_{sample}.ann.vcf.gz'),
	output:
		tsv = join(VCF, DATASET, OLD, 'variants', 'annotated', 'resistance', '{ref}_{drug}_{sample}.tsv'),
	conda:
		"../envs/alignment.yml"
	shell:
		"""
		bcftools filter -Oz -sFAIL -g3 -G10 -e '%QUAL<30 || FORMAT/DP<4' {input.vcf} > {params.tmp}
		bcftools index -f {params.tmp}
		bcftools query -R {params.bed}  -f '[%SAMPLE]\t%POS\t[%GT]\t%REF\t%ALT{{0}}\t%TYPE\t%QUAL\t%FILTER\t%INFO/DP\t[%INFO/DP4]\t[%INFO/ANN]\n' -i 'GT="alt"' {params.tmp} | \
			perl -p -e 's/\|/\t/g' | perl -p -e 's|{params.string}(.+).bam|$1|' > {output.tsv} 
		rm {params.tmp}
		"""

rule stats_coverage_old:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	params:
		ref	= REF
	input:
		expand(join(ALIGNMENTS, DATASET, OLD, REF + '_' + '{sample}.bam'), sample=SAMPLES)
	output:
		join(RESULTS, DATASET, OLD, 'stats.coverage.' + REF + '.txt'),
	conda:
		"../envs/alignment-0.1.20.yml"
	shell:
		"""
		for i in {input}; do echo -e $i","$(samtools depth -a $i | awk '{{sum+=$3}} END {{print sum/NR}}') | perl -pe 's/.+\/{params.ref}_(.+).bam/$1/'; done > {output}
		"""

rule stats_read_count_old:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	params:
		dir	= DATA + '\/' + DATASET +'\/' + OLD,
	input:
##		expand(join(DATA, DATASET, '{sample}_1.fastq.gz'), sample=SAMPLES)
		sampledata['forward'],
	output:
		join(RESULTS, DATASET, OLD, 'stats.read-count.txt'),
	conda:
		"../envs/alignment-0.1.20.yml"
	shell:
		"""
		for i in {input}; do echo -e $i","$(($(zcat $i | wc -l) / 4)) | perl -pe 's/{params.dir}\/(.+)_1.fastq.gz/$1/'; done > {output}
		"""

rule stats_combined_old:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	params:
		dir	= DATA + '\/' + DATASET + '\/' + OLD,
	input:
		reads = rules.stats_read_count.output,
		coverage = rules.stats_coverage.output,
	output:
		join(RESULTS, DATASET, OLD, 'stats.txt'),
	run:
		df_meta = pd.read_csv(meta_file, sep='\t')

		cols = ['file', 'counts']
		df_counts = pd.read_csv(str(input.reads), names=cols)
		
		cols = ['sample', 'coverage']
		df_coverage = pd.read_csv(str(input.coverage), names=cols)

		df_new = pd.merge(df_meta, df_counts, left_on='forward', right_on='file')
		df_all = pd.merge(df_new, df_coverage, on='sample')
			
		cols = ['sample', 'counts', 'coverage']
		df_all.to_csv(output[0], sep='\t', columns=cols, index=False)

rule lineages_old:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	input:
		expand(join(VCF, DATASET, OLD, REF + '_' + '{sample}.all.vcf.gz'), sample=SAMPLES)
	output:
		join(RESULTS, DATASET, OLD, 'lineages.txt'),
	conda:
		"alignment-bwa.yml"
	shell:
		"""
		for i in {input}; do perl scripts/check_lineages.pl scripts/lineages.txt $i; done > {output}
		"""

rule check_snps_all_old:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	params:
		dir	= DATA + '\/' + DATASET + '\/' + OLD,
	input:
		expand(join(VCF, DATASET, OLD, 'variants/annotated', REF + '_' + '{sample}.ann.vcf.gz'), sample=SAMPLES)
	output:
		join(RESULTS, DATASET, OLD, 'gene-snps-all.tsv'),
	conda:
		"alignment-bwa.yml"
	shell:
		"""
		for i in {input}; \
		do \
			j=`echo $i | perl -p -e 's/{params.dir}\/NC_000962_(.+).ann.vcf.gz/$1/'`; \
			echo "Sample: $j [File: $i]"; \
			echo $'\n'; \
			tabix -R scripts/snps.Chromosome-all.bed $i; \
			echo $'\n==============================\n'; \
		done > {output}
		"""

rule check_snps_bdq_old:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	params:
		dir	= DATA + '\/' + DATASET + '\/' + OLD,
	input:
		expand(join(VCF, DATASET, OLD, 'variants/annotated', REF + '_' + '{sample}.ann.vcf.gz'), sample=SAMPLES)
	output:
		join(RESULTS, DATASET, OLD, 'gene-snps-BDQ.txt'),
	conda:
		"alignment-bwa.yml"
	shell:
		"""
		for i in {input}; \
		do \
			j=`echo $i | perl -p -e 's/{params.dir}\/NC_000962_(.+).ann.vcf.gz/$1/'`; \
			echo "Sample: $j [File: $i]"; \
			echo $'\n'; \
			tabix -R scripts/snps.Chromosome-BDQ.bed $i; \
			echo $'\n==============================\n'; \
		done > {output}
		"""

rule check_snps_ptm_old:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	params:
		dir	= DATA + '\/' + DATASET + '\/' + OLD,
	input:
		expand(join(VCF, DATASET, OLD, 'variants/annotated', REF + '_' + '{sample}.ann.vcf.gz'), sample=SAMPLES)
	output:
		join(RESULTS, DATASET, OLD, 'gene-snps-PTM.txt'),
	conda:
		"alignment-bwa.yml"
	shell:
		"""
		for i in {input}; \
		do \
			j=`echo $i | perl -p -e 's/{params.dir}\/NC_000962_(.+).ann.vcf.gz/$1/'`; \
			echo "Sample: $j [File: $i]"; \
			echo $'\n'; \
			tabix -R scripts/snps.Chromosome-PTM.bed $i; \
			echo $'\n==============================\n'; \
		done > {output}
		"""

rule AMRPredict_old:
	threads:
		1
	params:
		dir	= DATA + '\/' + DATASET,
	input:
		expand(join(VCF, DATASET, REF + '_' + '{sample}.all.vcf.gz'), sample=SAMPLES)
	output:
		join(RESULTS, DATASET, 'AMRPredict.txt'),
	conda:
		"alignment-bwa.yml"
	shell:
		"""
		for i in {input}; \
		do \
			echo $i;	
			perl  /homedirs8/bunit/NGS/Projects/Mycobacterium/ResistDB/amrpredict.pl \
				--chrom NC_000962.3 \
				--gff /homedirs8/bunit/NGS/Projects/Mycobacterium/ResistDB/NC_000962.gff \
				--vcf $i; \
		done > {output}
		"""
