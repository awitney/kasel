rule versions:
	threads:
		1
	output:
		join(VERSIONS, 'alignment.txt')
	log:
		join(LOGS, 'versions.log')
	message:
		"Running alignment versions"
	conda:
		"../envs/alignment.yml"
	shell:
		"""
		( bwa &> {output}
		samtools --version &>> {output}  ) &> {log}
		"""

rule alignment_pe:
	threads:
		config['alignment_pe']['threads']
	resources:
		memory = config['alignment_pe']['memory']
	params:
		tmp = join(TMP, 'temp.{ref}_{sample}'),
	input:
		genome = join(GENOMES, '{ref}.fna'),
		r1 = lambda wildcards: get_seq(wildcards, 'forward'),
		r2 = lambda wildcards: get_seq(wildcards, 'reverse'),
	output:
		bam = join(ALIGNMENTS, '{ref}_{sample}.bam'),
		versions = join(ALIGNMENTS, '{ref}_{sample}.versions.txt'),
	log:
		join(LOGS, 'alignment_pe.{ref}.{sample}.log')
	message:
		"Running alignment_pe for {wildcards.sample} mapping to {wildcards.ref}"
	conda:
		"../envs/alignment.yml"
	shell:
		"""
		( bwa mem -t {threads} {input.genome} {input.r1} {input.r2} | samtools view -bS - | samtools sort -T {params.tmp} - | samtools rmdup - - > {output.bam} ) 2> {log}
		samtools index {output.bam}
		( echo "bwa " $(( bwa 2>&1 ) | grep Version ) > {output.versions} ) || exit 0
		samtools --version | grep "samtools " &>> {output.versions}
		"""

#		bwa &> {output.versions}
rule site_calling:
	threads:
		config['site_calling']['threads']
	resources:
		memory = config['site_calling']['memory']
	input:
		genome = join(GENOMES, '{ref}.fna'),
		bam = join(ALIGNMENTS, '{ref}_{sample}.bam'),
	output:
		vcf = join(VCF, '{ref}_{sample}.all.vcf.gz'),
		versions = join(VCF, '{ref}_{sample}.all.vcf.versions.txt'),
	log:
		join(LOGS, 'site_calling.{ref}.{sample}.log')
	message:
		"Running site_calling for {wildcards.sample} mapped to {wildcards.ref}"
	conda:
		"../envs/alignment.yml"
	shell:
		"""
		( bcftools mpileup -Ou -f {input.genome} {input.bam} --annotate 'FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR' | \
			bcftools call --threads 1 -c | \
			bcftools filter -sMixed -e '(DP4[0]+DP4[1])>10 & (DP4[2]+DP4[3])>10' - | \
			bcftools filter -sDepth -e 'FORMAT/DP<10' - | \
			bcftools filter -sLowQual -g3 -G10 -e '%QUAL<30 || RPB<0.1' - | \
			bgzip > {output.vcf} && tabix -p vcf {output.vcf} ) 2> {log}
		bcftools --version | grep "bcftools " &>> {output.versions}
		"""

rule variant_calling:
	threads:
		config['variant_calling']['threads']
	resources:
		memory = config['variant_calling']['memory']
	params:
		ref = '{ref}.3'
	input:
		genome = join(GENOMES, '{ref}.fna'),
		bam = join(ALIGNMENTS, '{ref}_{sample}.bam'),
	output:
		vcf = join(VCF, 'variants', '{ref}_{sample}.vcf.gz'),
		versions = join(VCF, '{ref}_{sample}.vcf.versions.txt'),
	log:
		join(LOGS, 'variant_calling.{ref}.{sample}.log')
	message:
		"Running variant_calling for {wildcards.sample} mapped to {wildcards.ref}"
	conda:
		"../envs/alignment.yml"
	shell:
		"""
		( bcftools mpileup -Ou -f {input.genome} {input.bam} --annotate 'FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR' | \
			bcftools call --threads 1 -cv | \
			bcftools filter -sMixed -e '(DP4[2]+DP4[3])/sum(DP4)<0.75 | (DP4[0]+DP4[1])>10' - | \
			bcftools filter -sDepth -e 'FORMAT/DP<10' - | \
			bcftools filter -sLowQual -g3 -G10 -e '%QUAL<30 || RPB<0.1' - | \
			perl -p -e 's/^{params.ref}/Chromosome/' | bgzip > {output.vcf} && tabix -p vcf {output.vcf} ) 2> {log}
		bcftools --version | grep "bcftools " &> {output.versions}
		"""

rule snp_annotation:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	input:
		join(VCF, 'variants', '{ref}_{sample}.vcf.gz'),
	output:
		vcf = join(VCF, 'variants', 'annotated', '{ref}_{sample}.ann.vcf.gz'),
		versions = join(VCF, 'variants', 'annotated', '{ref}_{sample}.ann.vcf.versions.txt'),
	log:
		join(LOGS, 'snp_annotation.{ref}.{sample}.log')
	message:
		"Running snp_annotation for {wildcards.sample} mapped to {wildcards.ref}"
	conda:
		"../envs/alignment.yml"
	shell:
		"""
        ( snpEff eff -no-downstream -no-upstream -no-utr -o vcf Mycobacterium_tuberculosis_h37rv {input} | bgzip > {output.vcf} && bcftools index {output.vcf} ) &> {log}
		snpEff -version &> {output.versions} 
		"""

rule snp_report_all:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	input:
		tsv = expand(join(VCF, 'variants/annotated', '{ref}_{sample}.tsv'), ref=REF, sample=SAMPLES),
	output:
		tsv = join(RESULTS, 'variants', '{ref}_' + DATASET + '.tsv'),
	message:
		"Running snp_report_all for {wildcards.ref}"
	conda:
		"../envs/alignment.yml"
	shell:
		"""
		echo "Sample\tPos\tGT\tREF\tALT\tType\tQUAL\tFilter\tDP\tDP4\tANN\t" > {output.tsv}
		cat {input.tsv} >> {output.tsv}
		"""

rule snp_report:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	params:
		string = join(ALIGNMENTS, '{ref}_')
	input:
		vcf = join(VCF, 'variants/annotated', '{ref}_{sample}.ann.vcf.gz'),
	output:
		tsv = join(VCF, 'variants/annotated', '{ref}_{sample}.tsv'),
		versions = join(VCF, 'variants/annotated', '{ref}_{sample}.versions.txt'),
	log:
		join(LOGS, 'snp_report.{ref}.{sample}.log')
	message:
		"Running snp_report for {wildcards.sample} mapped to {wildcards.ref}"
	conda:
		"../envs/alignment.yml"
	shell:
		"""
		( bcftools index -f {input.vcf}
		bcftools filter -sFAIL -g3 -G10 -e '%QUAL<30 || FORMAT/DP<4' {input.vcf} | \
			bcftools query -f '[%SAMPLE]\t%POS\t[%GT]\t%REF\t%ALT{{0}}\t%TYPE\t%QUAL\t%FILTER\t%INFO/DP\t[%INFO/DP4]\t[%INFO/ANN]\n' -i 'GT="alt"' - | \
			perl -p -e 's/\|/\t/g' > {output.tsv} ) &> {log}
		bcftools --version | grep "bcftools " &> {output.versions}
		"""
#			perl -p -e 's/\|/\t/g' | perl -p -e 's|{params.string}(.+).bam|$1|' > {output.tsv} 

rule snp_report_resistance_all:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	input:
		tsv = expand(join(VCF, 'variants/annotated/resistance', '{ref}_{{drug}}_{sample}.tsv'), ref=REF, sample=SAMPLES),
	output:
		tsv = join(RESULTS, 'variants/resistance', '{ref}_' + DATASET + '_{drug}.tsv'),
	message:
		"Running snp_report_resistance_all for {wildcards.ref} drug {wildcards.drug}"
	conda:
		"../envs/alignment.yml"
	shell:
		"""
		echo "Sample\tPos\tGT\tREF\tALT\tType\tQUAL\tFilter\tDP\tDP4\tANN\t" > {output.tsv}
		cat {input.tsv} >> {output.tsv}
		"""

rule snp_report_resistance_all_summary:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	params:
		finddir = join(VCF, 'variants/annotated/resistance/'),
		string = join(VCF, 'variants/annotated/resistance/' + REF + '_'),
	input:
		tsv1 = expand(join(RESULTS, 'variants/resistance', '{ref}_' + DATASET + '_BDQ.tsv'), ref=REF),
		tsv2 = expand(join(RESULTS, 'variants/resistance', '{ref}_' + DATASET + '_PTM.tsv'), ref=REF),
		tsv3 = expand(join(RESULTS, 'variants/resistance', '{ref}_' + DATASET + '_LZD.tsv'), ref=REF),
	output:
		out = join(RESULTS, 'variants/resistance', '{ref}_' + DATASET + '.tsv'),
	message:
		"Running snp_report_resistance_all_summary for {wildcards.ref}"
	conda:
		"../envs/alignment.yml"
	shell:
		"""
		wc -l  $(find {params.finddir} -name "*_BDQ*.tsv") | perl -p -e 's:[ ]+(\d+)[ ]+{params.string}(PTM|BDQ|LZD)_(.+).tsv:$3\t$2\t$1:' | sort -k2 > {output.out}
		wc -l  $(find {params.finddir} -name "*_PTM*.tsv") | perl -p -e 's:[ ]+(\d+)[ ]+{params.string}(PTM|BDQ|LZD)_(.+).tsv:$3\t$2\t$1:' | sort -k2 >> {output.out}
		wc -l  $(find {params.finddir} -name "*_LZD*.tsv") | perl -p -e 's:[ ]+(\d+)[ ]+{params.string}(PTM|BDQ|LZD)_(.+).tsv:$3\t$2\t$1:' | sort -k2 >> {output.out}
		"""

rule snp_annotation_filter:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	input:
		vcf = join(VCF, 'variants', 'annotated', '{ref}_{sample}.ann.vcf.gz'),
	output:
		tmp = join(VCF, 'variants', 'annotated', '{ref}_{sample}.ann.vcf.tmp.gz'),
		csi = join(VCF, 'variants', 'annotated', '{ref}_{sample}.ann.vcf.tmp.gz.csi'),
		versions = join(VCF, 'variants', 'annotated', '{ref}_{sample}.ann.vcf.tmp.versions.txt'),
	log:
		join(LOGS, 'snp_annotation_filter.{ref}.{sample}.log')
	message:
		"Running snp_annotation_filter for {wildcards.sample} mapped to {wildcards.ref}"
	conda:
		"../envs/alignment.yml"
	shell:
		"""
		echo "Starting snp_annotation_filter" > {log}
		bcftools filter -Oz -sFAIL -g3 -G10 -e '%QUAL<30 || FORMAT/DP<4' {input.vcf} > {output.tmp} 2>> {log};
		bcftools index -f {output.tmp} 2>> {log};
		echo "Finished snp_annotation_filter" >> {log}
		bcftools --version | grep "bcftools " &> {output.versions}
		"""

rule snp_report_resistance:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	params:
		string = join(ALIGNMENTS, '{ref}_'),
		bedfile = join(config['kasel-data'], 'snps.Chromosome-{drug}.bed'),
	input:
		tmp = join(VCF, 'variants', 'annotated', '{ref}_{sample}.ann.vcf.tmp.gz'),
		csi = join(VCF, 'variants', 'annotated', '{ref}_{sample}.ann.vcf.tmp.gz.csi'),
	output:
		tsv = join(VCF, 'variants', 'annotated', 'resistance', '{ref}_{drug}_{sample}.tsv'),
		versions = join(VCF, 'variants', 'annotated', 'resistance', '{ref}_{drug}_{sample}.versions.txt'),
	log:
		join(LOGS, 'snp_report_resistance.{ref}.{drug}.{sample}.log')
	message:
		"Running snp_report_resistance for {wildcards.sample} mapped to {wildcards.ref} drug {wildcards.drug}"
	conda:
		"../envs/alignment.yml"
	shell:
		"""
		echo "Starting snp_report_resistance" > {log}
		bcftools query -R {params.bedfile}  -f '[%SAMPLE]\\t%POS\\t[%GT]\\t%REF\\t%ALT{{0}}\\t%TYPE\\t%QUAL\\t%FILTER\\t%INFO/DP\\t[%INFO/DP4]\\t[%INFO/ANN]\\n' -i 'GT="alt"' {input.tmp} | \
			perl -p -e 's/\|/\t/g' | perl -p -e 's|{params.string}(.+).bam|$1|' > {output.tsv} 2>> {log};
		echo "Finished snp_report_resistance" >> {log}
		bcftools --version | grep "bcftools " &> {output.versions}
		"""

rule stats_coverage:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	params:
		ref	= REF
	input:
		expand(join(ALIGNMENTS, REF + '_' + '{sample}.bam'), sample=SAMPLES)
	output:
		join(RESULTS, 'stats.coverage.' + REF + '.txt'),
	log:
		join(LOGS, 'stats_coverage.log')
	message:
		"Running stats_coverage"
	conda:
		"../envs/alignment.yml"
	shell:
		"""
		echo "Starting stats_coverage" > {log}
		for i in {input}; do echo -e $i","$(samtools depth -a $i | awk '{{sum+=$3}} END {{print sum/NR}}') | perl -pe 's/.+\/{params.ref}_(.+).bam/$1/'; done > {output} 2>> {log}
		echo "Finishing stats_coverage" >> {log}
		"""

rule stats_read_count:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	params:
		dir	= DATA,
	input:
		sampledata['forward'],
	output:
		join(RESULTS, 'stats.read-count.txt'),
	log:
		join(LOGS, 'stats_read_count.log')
	message:
		"Running stats_read_count"
	conda:
		"../envs/alignment.yml"
	shell:
		"""
		echo "Starting stats_read_count" > {log}
		for i in {input}; do echo -e $i","$(($(zcat $i | wc -l) / 4)) | perl -pe 's/.+\/{params.dir}\/(.+)_1.fastq.gz/$1/'; done > {output} 2>> {log}
		echo "Finishing stats_read_count" >> {log}
		"""

rule stats_combined:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	params:
		dir	= DATA,
	input:
		reads = rules.stats_read_count.output,
		coverage = rules.stats_coverage.output,
	output:
		join(RESULTS, 'stats.txt'),
	message:
		"Running stats_combined"
	run:
		df_samples = pd.read_csv(samples_file, sep='\t')

		cols = ['file', 'counts']
		df_counts = pd.read_csv(str(input.reads), names=cols)
		
		cols = ['sample', 'coverage']
		df_coverage = pd.read_csv(str(input.coverage), names=cols)

		df_new = pd.merge(df_samples, df_counts, left_on='forward', right_on='file')
		df_all = pd.merge(df_new, df_coverage, on='sample')
			
		cols = ['sample', 'counts', 'coverage']
		df_all.to_csv(output[0], sep='\t', columns=cols, index=False)

rule lineages:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	input:
		expand(join(VCF, REF + '_' + '{sample}.all.vcf.gz'), sample=SAMPLES)
	output:
		join(RESULTS, 'lineages.txt'),
	message:
		"Running lineages"
	conda:
		"alignment-bwa.yml"
	shell:
		"""
		for i in {input}; do perl scripts/check_lineages.pl scripts/lineages.txt $i; done > {output}
		"""

rule check_snps_all:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	params:
		dir	= DATA,
		bedfile = join(config['kasel-data'], 'snps.Chromosome-all.bed'),
	input:
		expand(join(VCF, 'variants/annotated', REF + '_' + '{sample}.ann.vcf.gz'), sample=SAMPLES)
	output:
		join(RESULTS, 'gene-snps-all.tsv'),
	log:
		join(LOGS, 'check_snps_all.log')
	message:
		"Running check_snps_all"
	conda:
		"../envs/alignment.yml"
	shell:
		"""
		for i in {input}; \
		do \
			j=`echo $i | perl -p -e 's/{params.dir}\/NC_000962_(.+).ann.vcf.gz/$1/'`; \
			echo "Sample: $j [File: $i]"; \
			echo $'\n'; \
			tabix -R {params.bedfile} $i; \
			echo $'\n==============================\n'; \
		done > {output} 2> {log}
		"""

rule check_snps_bdq:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	params:
		dir	= DATA,
		bedfile = join(config['kasel-data'], 'snps.Chromosome-BDQ.bed'),
	input:
		samples = expand(join(VCF, 'variants/annotated', REF + '_' + '{sample}.ann.vcf.gz'), sample=SAMPLES),
	output:
		join(RESULTS, 'gene-snps-BDQ.txt'),
	log:
		join(LOGS, 'check_snps_bdq.log')
	message:
		"Running check_snps_bdq"
	conda:
		"../envs/alignment.yml"
	shell:
		"""
		echo {params.bedfile} > {log};
		less {params.bedfile} >> {log};
		for i in {input.samples}; \
		do \
			j=`echo $i | perl -p -e 's/{params.dir}\/NC_000962_(.+).ann.vcf.gz/$1/'`; \
			echo "Sample: $j [File: $i]"; \
			echo $'\n'; \
			tabix -R {params.bedfile} $i; \
			echo $'\n==============================\n'; \
		done > {output} 2>> {log}
		"""

rule check_snps_ptm:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	params:
		dir	= DATA,
		bedfile = join(config['kasel-data'], 'snps.Chromosome-PTM.bed'),
	input:
		expand(join(VCF, 'variants/annotated', REF + '_' + '{sample}.ann.vcf.gz'), sample=SAMPLES)
	output:
		join(RESULTS, 'gene-snps-PTM.txt'),
	log:
		join(LOGS, 'check_snps_ptm.log')
	message:
		"Running check_snps_ptm"
	conda:
		"../envs/alignment.yml"
	shell:
		"""
		echo {params.bedfile} > {log};
		less {params.bedfile} >> {log};
		for i in {input}; \
		do \
			j=`echo $i | perl -p -e 's/{params.dir}\/NC_000962_(.+).ann.vcf.gz/$1/'`; \
			echo "Sample: $j [File: $i]"; \
			echo $'\n'; \
			tabix -R {params.bedfile} $i; \
			echo $'\n==============================\n'; \
		done > {output}
		"""

rule check_snps_lzd:
	threads:
		config['default']['threads']
	resources:
		memory = config['default']['memory']
	params:
		dir	= DATA,
		bedfile = join(config['kasel-data'], 'snps.Chromosome-LZD.bed'),
	input:
		expand(join(VCF, 'variants/annotated', REF + '_' + '{sample}.ann.vcf.gz'), sample=SAMPLES)
	output:
		join(RESULTS, 'gene-snps-LZD.txt'),
	log:
		join(LOGS, 'check_snps_lzd.log')
	message:
		"Running check_snps_lzd"
	conda:
		"../envs/alignment.yml"
	shell:
		"""
		echo {params.bedfile} > {log};
		less {params.bedfile} >> {log};
		for i in {input}; \
		do \
			j=`echo $i | perl -p -e 's/{params.dir}\/NC_000962_(.+).ann.vcf.gz/$1/'`; \
			echo "Sample: $j [File: $i]"; \
			echo $'\n'; \
			tabix -R {params.bedfile} $i; \
			echo $'\n==============================\n'; \
		done > {output}
		"""

rule AMRPredict:
	threads:
		1
	params:
		dir	= DATA,
	input:
		expand(join(VCF, REF + '_' + '{sample}.all.vcf.gz'), sample=SAMPLES)
	output:
		join(RESULTS, 'AMRPredict.txt'),
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
