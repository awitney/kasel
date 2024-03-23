
rule alignment_pe_legacy:
	threads:
		config['alignment_pe']['threads']
	resources:
		memory = config['alignment_pe']['memory']
	params:
		tmp = join(TMP, LEGACY, 'temp.{ref}_{sample}'),
	input:
		genome = join(GENOMES, '{ref}.fna'),
		r1 = lambda wildcards: get_seq(wildcards, 'forward'),
		r2 = lambda wildcards: get_seq(wildcards, 'reverse'),
	output:
		bam = join(ALIGNMENTS, LEGACY, '{ref}_{sample}.bam'),
		versions = join(ALIGNMENTS, LEGACY, '{ref}_{sample}.versions.txt'),
	log:
		join(LOGS, LEGACY, 'alignment_pe.{ref}.{sample}.log')
	message:
		"Running alignment_pe_legacy for {wildcards.sample} mapping to {wildcards.ref}"
	conda:
		"../envs/alignment-legacy.yml"
	shell:
		"""
		( echo "bwa " $(( bwa 2>&1 ) | grep Version ) > {output.versions} ) || exit 0
		( echo "samtools " $(( samtools 2>&1 ) | grep Version ) >> {output.versions} ) || exit 0
		( bwa mem -t {threads} {input.genome} {input.r1} {input.r2} | samtools view -bS - | samtools sort -o - {params.tmp} - | samtools rmdup - - > {output.bam} ) 2> {log}
		samtools index {output.bam}
		"""

rule site_calling_legacy:
	threads:
		config['site_calling']['threads']
	resources:
		memory = config['site_calling']['memory']
	input:
		genome = join(GENOMES, '{ref}.fna'),
		bam = join(ALIGNMENTS, LEGACY, '{ref}_{sample}.bam'),
	output:
		vcf = join(VCF, LEGACY, '{ref}_{sample}.all.vcf.gz'),
		versions = join(VCF, LEGACY, '{ref}_{sample}.all.vcf.versions.txt'),
	log:
		join(LOGS, LEGACY, 'site_calling.{ref}.{sample}.log')
	message:
		"Running site_calling_legacy for {wildcards.sample} mapped to {wildcards.ref}"
	conda:
		"../envs/alignment-legacy.yml"
	shell:
		"""
		( echo "samtools " $(( samtools 2>&1 ) | grep Version ) >> {output.versions} ) || exit 0
		( samtools mpileup -gf {input.genome} {input.bam} | bcftools view -cg - | bgzip > {output.vcf} && tabix -p vcf {output.vcf} ) &> {log}
		"""

rule variant_calling_legacy:
	threads:
		config['variant_calling']['threads']
	resources:
		memory = config['variant_calling']['memory']
	params:
		ref = '{ref}.3',
		bcf = join(VCF, LEGACY, '{ref}_{sample}.bcf')
	input:
		genome = join(GENOMES, '{ref}.fna'),
		bam = join(ALIGNMENTS, LEGACY, '{ref}_{sample}.bam'),
	output:
		vcf = join(VCF, LEGACY, 'variants', '{ref}_{sample}.vcf.gz'),
		versions = join(VCF, LEGACY, '{ref}_{sample}.vcf.versions.txt'),
	log:
		join(LOGS, LEGACY, 'variant_calling.{ref}.{sample}.log')
	message:
		"Running variant_calling_legacy for {wildcards.sample} mapped to {wildcards.ref}"
	conda:
		"../envs/alignment-legacy.yml"
	shell:
		"""
		( echo "samtools " $(( samtools 2>&1 ) | grep Version ) >> {output.versions} ) || exit 0
		( samtools mpileup -ugf {input.genome} {input.bam} | bcftools view -bvcg - > {params.bcf} && bcftools view {params.bcf} | vcfutils.pl varFilter -D 20000 | sed 's/NC_000962.3/Chromosome/g' | bgzip > {output.vcf} && tabix -p vcf {output.vcf} ) &> {log}
		rm {params.bcf}
		"""

rule snp_annotation_legacy:
	threads:
		config['snp_annotation_legacy']['threads']
	resources:
		memory = config['snp_annotation_legacy']['memory']
	input:
		join(VCF, LEGACY, 'variants', '{ref}_{sample}.vcf.gz'),
	output:
		vcf = join(VCF, LEGACY, 'variants', 'annotated', '{ref}_{sample}.ann.vcf.gz'),
		versions = join(VCF, LEGACY, 'variants', 'annotated', '{ref}_{sample}.ann.vcf.versions.txt'),
	log:
		join(LOGS, LEGACY, 'snp_annotation_legacy.{ref}.{sample}.log')
	message:
		"Running snp_annotation_legacy for {wildcards.sample} mapped to {wildcards.ref}"
	conda:
		"../envs/alignment-legacy.yml"
	shell:
		"""
		snpEff -version &> {output.versions} 
        ( snpEff eff -no-downstream -no-upstream -no-utr -o vcf Mycobacterium_tuberculosis_h37rv {input} | bgzip > {output.vcf} && bcftools index {output.vcf} ) &> {log}
		"""

rule snp_report_all_legacy:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	input:
		tsv = expand(join(VCF, LEGACY, 'variants/annotated', '{ref}_{sample}.tsv'), ref=REF, sample=SAMPLES),
	output:
		tsv = join(RESULTS, LEGACY, 'variants', '{ref}_' + DATASET + '.tsv'),
	message:
		"Running snp_report_all_legacy for {wildcards.ref}"
	conda:
		"../envs/alignment-legacy.yml"
	shell:
		"""
		echo "Sample\tPos\tGT\tREF\tALT\tType\tQUAL\tFilter\tDP\tDP4\tANN\t" > {output.tsv}
		perl -p -e 's/.+\/NC_000962_(.+).bam/$1/' {input.tsv} >> {output.tsv}
		"""

rule snp_report_legacy:
	threads:
		config['snp_report_legacy']['threads']
	resources:
		memory = config['snp_report_legacy']['memory']
	params:
		string = join(ALIGNMENTS, LEGACY, '{ref}_')
	input:
		vcf = join(VCF, LEGACY, 'variants/annotated', '{ref}_{sample}.ann.vcf.gz'),
	output:
		tsv = join(VCF, LEGACY, 'variants/annotated', '{ref}_{sample}.tsv'),
		versions = join(VCF, LEGACY, 'variants/annotated', '{ref}_{sample}.versions.txt'),
	log:
		join(LOGS, LEGACY, 'snp_report.{ref}.{sample}.log')
	message:
		"Running snp_report_legacy for {wildcards.sample} mapped to {wildcards.ref}"
	conda:
		"../envs/alignment.yml"
	shell:
		"""
		bcftools --version | grep "bcftools " &> {output.versions}
		( bcftools index -f {input.vcf}
		bcftools filter -sFAIL -g3 -G10 -e '%QUAL<30 || FORMAT/DP<4' {input.vcf} | \
			bcftools query -f '[%SAMPLE]\\t%POS\\t[%GT]\\t%REF\\t%ALT{{0}}\\t%TYPE\\t%QUAL\\t%FILTER\\t%INFO/DP\\t[%INFO/DP4]\\t[%INFO/ANN]\\n' -i 'GT="alt"' - | \
			perl -p -e 's/\|/\t/g' > {output.tsv} ) &> {log}
		"""

rule snp_report_resistance_all_legacy:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	input:
		tsv = expand(join(VCF, LEGACY, 'variants/annotated/resistance', '{ref}_{{drug}}_{sample}.tsv'), ref=REF, sample=SAMPLES),
	output:
		tsv = join(RESULTS, LEGACY, 'variants/resistance', '{ref}_' + DATASET + '_{drug}.tsv'),
	message:
		"Running snp_report_resistance_all_legacy for {wildcards.ref} drug {wildcards.drug}"
	conda:
		"../envs/alignment.yml"
	shell:
		"""
		echo "Sample\tPos\tGT\tREF\tALT\tType\tQUAL\tFilter\tDP\tDP4\tANN\t" > {output.tsv}
		cat {input.tsv} >> {output.tsv}
		"""

rule snp_report_resistance_all_summary_legacy:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	params:
		finddir = join(VCF, LEGACY, 'variants/annotated/resistance/'),
		string = join(VCF, LEGACY, 'variants/annotated/resistance/' + REF + '_'),
	input:
		tsv1 = expand(join(RESULTS, LEGACY, 'variants/resistance', '{ref}_' + DATASET + '_BDQ.tsv'), ref=REF),
		tsv2 = expand(join(RESULTS, LEGACY, 'variants/resistance', '{ref}_' + DATASET + '_PTM.tsv'), ref=REF),
		tsv3 = expand(join(RESULTS, LEGACY, 'variants/resistance', '{ref}_' + DATASET + '_LZD.tsv'), ref=REF),
	output:
		out = join(RESULTS, LEGACY, 'variants/resistance', '{ref}_' + DATASET + '.tsv'),
	message:
		"Running snp_report_resistance_all_summary_legacy for {wildcards.ref}"
	conda:
		"../envs/alignment.yml"
	shell:
		"""
		wc -l  $(find {params.finddir} -name "*_BDQ*.tsv") | perl -p -e 's:[ ]+(\d+)[ ]+{params.string}(PTM|BDQ|LZD)_(.+).tsv:$3\t$2\t$1:' | sort -k2 > {output.out}
		wc -l  $(find {params.finddir} -name "*_PTM*.tsv") | perl -p -e 's:[ ]+(\d+)[ ]+{params.string}(PTM|BDQ|LZD)_(.+).tsv:$3\t$2\t$1:' | sort -k2 >> {output.out}
		wc -l  $(find {params.finddir} -name "*_LZD*.tsv") | perl -p -e 's:[ ]+(\d+)[ ]+{params.string}(PTM|BDQ|LZD)_(.+).tsv:$3\t$2\t$1:' | sort -k2 >> {output.out}
		"""

rule snp_annotation_filter_legacy:
	threads:
		config['snp_annotation_filter_legacy']['threads']
	resources:
		memory = config['snp_annotation_filter_legacy']['memory']
	input:
		vcf = join(VCF, LEGACY, 'variants', 'annotated', '{ref}_{sample}.ann.vcf.gz'),
	output:
		tmp = join(VCF, LEGACY, 'variants', 'annotated', '{ref}_{sample}.ann.vcf.tmp.gz'),
		csi = join(VCF, LEGACY, 'variants', 'annotated', '{ref}_{sample}.ann.vcf.tmp.gz.csi'),
		versions = join(VCF, LEGACY, 'variants', 'annotated', '{ref}_{sample}.ann.vcf.tmp.versions.txt'),
	log:
		join(LOGS, 'snp_annotation_filter_legacy.{ref}.{sample}.log')
	message:
		"Running snp_annotation_filter_legacy for {wildcards.sample} mapped to {wildcards.ref}"
	conda:
		"../envs/alignment.yml"
	shell:
		"""
		bcftools --version | grep "bcftools " &> {output.versions}
		echo "Starting snp_annotation_filter_legacy" > {log}
		bcftools filter -Oz -sFAIL -g3 -G10 -e '%QUAL<30 || FORMAT/DP<4' {input.vcf} > {output.tmp} 2>> {log};
		bcftools index -f {output.tmp} 2>> {log};
		echo "Finished snp_annotation_filter_legacy" >> {log}
		"""

rule snp_report_resistance_legacy:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	params:
		string = join(ALIGNMENTS, LEGACY, '{ref}_'),
		bedfile = join(config['kasel-data'], 'snps.Chromosome-{drug}.bed'),
	input:
		tmp = join(VCF, LEGACY, 'variants', 'annotated', '{ref}_{sample}.ann.vcf.tmp.gz'),
		csi = join(VCF, LEGACY, 'variants', 'annotated', '{ref}_{sample}.ann.vcf.tmp.gz.csi'),
	output:
		tsv = join(VCF, LEGACY, 'variants', 'annotated', 'resistance', '{ref}_{drug}_{sample}.tsv'),
		versions = join(VCF, LEGACY, 'variants', 'annotated', 'resistance', '{ref}_{drug}_{sample}.versions.txt'),
	log:
		join(LOGS, LEGACY, 'snp_report_resistance.{ref}.{drug}.{sample}.log')
	message:
		"Running snp_report_resistance_legacy for {wildcards.sample} mapped to {wildcards.ref} drug {wildcards.drug}"
	conda:
		"../envs/alignment.yml"
	shell:
		"""
		bcftools --version | grep "bcftools " &> {output.versions}
		echo "Starting snp_report_resistance_legacy" > {log}
		bcftools query -R {params.bedfile}  -f '[%SAMPLE]\t%POS\t[%GT]\t%REF\t%ALT{{0}}\t%TYPE\t%QUAL\t%FILTER\t%INFO/DP\t[%INFO/DP4]\t[%INFO/ANN]\n' -i 'GT="alt"' {input.tmp} | \
			perl -p -e 's/\|/\t/g' | perl -p -e 's|{params.string}(.+).bam|$1|' > {output.tsv} 2>> {log}
		echo "Finished snp_report_resistance_legacy" >> {log}
		"""

rule stats_coverage_legacy:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	params:
		ref	= REF
	input:
		expand(join(ALIGNMENTS, LEGACY, REF + '_' + '{sample}.bam'), sample=SAMPLES)
	output:
		join(RESULTS, LEGACY, 'stats.coverage.' + REF + '.txt'),
	message:
		"Running stats_coverage_legacy"
	conda:
		"../envs/alignment-legacy.yml"
	shell:
		"""
		for i in {input}; do echo -e $i","$(samtools depth -a $i | awk '{{sum+=$3}} END {{print sum/NR}}') | perl -pe 's/.+\/{params.ref}_(.+).bam/$1/'; done > {output}
		"""

rule stats_read_count_legacy:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	params:
		dir	= DATA + '\/' + LEGACY,
	input:
		sampledata['forward'],
	output:
		join(RESULTS, LEGACY, 'stats.read-count.txt'),
	message:
		"Running stats_read_count_legacy"
	conda:
		"../envs/alignment-legacy.yml"
	shell:
		"""
		for i in {input}; do echo -e $i","$(($(zcat $i | wc -l) / 4)) | perl -pe 's/{params.dir}\/(.+)_1.fastq.gz/$1/'; done > {output}
		"""

rule stats_combined_legacy:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	params:
		dir	= DATA + '\/' + LEGACY,
	input:
		reads = rules.stats_read_count.output,
		coverage = rules.stats_coverage.output,
	output:
		join(RESULTS, LEGACY, 'stats.txt'),
	message:
		"Running stats_combined_legacy"
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

rule lineages_legacy:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	input:
		expand(join(VCF, LEGACY, REF + '_' + '{sample}.all.vcf.gz'), sample=SAMPLES)
	output:
		join(RESULTS, LEGACY, 'lineages.txt'),
	message:
		"Running lineages_legacy"
	conda:
		"../envs/alignment-legacy.yml"
	shell:
		"""
		for i in {input}; do perl scripts/check_lineages.pl scripts/lineages.txt $i; done > {output}
		"""

rule check_snps_all_legacy:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	params:
		dir	= DATA + '\/' + LEGACY,
		bedfile = join(config['kasel-data'], 'snps.Chromosome-all.bed'),
	input:
		expand(join(VCF, LEGACY, 'variants/annotated', REF + '_' + '{sample}.ann.vcf.gz'), sample=SAMPLES)
	output:
		join(RESULTS, LEGACY, 'gene-snps-all.tsv'),
	message:
		"Running check_snps_all_legacy"
	conda:
		"../envs/alignment-legacy.yml"
	shell:
		"""
		for i in {input}; \
		do \
			j=`echo $i | perl -p -e 's/{params.dir}\/NC_000962_(.+).ann.vcf.gz/$1/'`; \
			echo "Sample: $j [File: $i]"; \
			echo $'\n'; \
			tabix -R {params.bedfile} $i; \
			echo $'\n==============================\n'; \
		done > {output}
		"""

rule check_snps_bdq_legacy:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	params:
		dir	= DATA + '\/' + LEGACY,
		bedfile = join(config['kasel-data'], 'snps.Chromosome-BDQ.bed'),
	input:
		expand(join(VCF, LEGACY, 'variants/annotated', REF + '_' + '{sample}.ann.vcf.gz'), sample=SAMPLES)
	output:
		join(RESULTS, LEGACY, 'gene-snps-BDQ.txt'),
	log:
		join(LOGS, LEGACY, 'check_snps_bdq.log')
	message:
		"Running check_snps_bdq_legacy"
	conda:
		"../envs/alignment-legacy.yml"
	shell:
		"""
		for i in {input}; \
		do \
			j=`echo $i | perl -p -e 's/{params.dir}\/NC_000962_(.+).ann.vcf.gz/$1/'`; \
			echo "Sample: $j [File: $i]"; \
			echo $'\n'; \
			tabix -R {params.bedfile} $i; \
			echo $'\n==============================\n'; \
		done > {output} 2>> {log}
		"""

rule check_snps_ptm_legacy:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	params:
		dir	= DATA + '\/' + LEGACY,
		bedfile = join(config['kasel-data'], 'snps.Chromosome-PTM.bed'),
	input:
		expand(join(VCF, LEGACY, 'variants/annotated', REF + '_' + '{sample}.ann.vcf.gz'), sample=SAMPLES)
	output:
		join(RESULTS, LEGACY, 'gene-snps-PTM.txt'),
	log:
		join(LOGS, LEGACY, 'check_snps_ptm.log')
	message:
		"Running check_snps_ptm"
	conda:
		"../envs/alignment-legacy.yml"
	shell:
		"""
		for i in {input}; \
		do \
			j=`echo $i | perl -p -e 's/{params.dir}\/NC_000962_(.+).ann.vcf.gz/$1/'`; \
			echo "Sample: $j [File: $i]"; \
			echo $'\n'; \
			tabix -R {params.bedfile} $i; \
			echo $'\n==============================\n'; \
		done > {output} 2>> {log}
		"""

rule check_snps_lzd_legacy:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	params:
		dir	= DATA + '\/' + LEGACY,
		bedfile = join(config['kasel-data'], 'snps.Chromosome-LZD.bed'),
	input:
		expand(join(VCF, LEGACY, 'variants/annotated', REF + '_' + '{sample}.ann.vcf.gz'), sample=SAMPLES)
	output:
		join(RESULTS, LEGACY, 'gene-snps-LZD.txt'),
	log:
		join(LOGS, LEGACY, 'check_snps_lzd.log')
	message:
		"Running check_snps_lzd"
	conda:
		"../envs/alignment-legacy.yml"
	shell:
		"""
		for i in {input}; \
		do \
			j=`echo $i | perl -p -e 's/{params.dir}\/NC_000962_(.+).ann.vcf.gz/$1/'`; \
			echo "Sample: $j [File: $i]"; \
			echo $'\n'; \
			tabix -R {params.bedfile} $i; \
			echo $'\n==============================\n'; \
		done > {output} 2>> {log}
		"""

rule AMRPredict_legacy:
	threads:
		1
	params:
		dir	= DATA,
	input:
		expand(join(VCF, REF + '_' + '{sample}.all.vcf.gz'), sample=SAMPLES)
	output:
		join(RESULTS, 'AMRPredict.txt'),
	conda:
		"../envs/alignment-legacy.yml"
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
