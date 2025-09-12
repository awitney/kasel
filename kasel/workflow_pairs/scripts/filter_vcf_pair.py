
import sys
import re
import argparse
import vcf

def get_args():

	parser = argparse.ArgumentParser(description='Process some integers.')

	parser.add_argument('--vcf1', type=str, required=True, help='VCF 1 filename')
	parser.add_argument('--vcf2', type=str, required=True, help='VCF 2 filename')
	parser.add_argument('--site', type=int, help='Show data for individual site')

	return parser

def main():

	args = get_args().parse_args()

	print("","Sample Details", "No Reads", "Coverage", sep="\t", file=sys.stderr)
	print("BASELINE (or equilivalent)", args.vcf1, sep="\t", file=sys.stderr)
	print("FOLLOW-UP", args.vcf2, sep="\t", file=sys.stderr)
	print("", sep="\t", file=sys.stderr)
	print("Phylogenetic No. of SNPs Different:", sep="\t", file=sys.stderr)
	print("", sep="\t", file=sys.stderr)
#	print("If SNPs different is <20 the table below listing the details of the SNPs must be complete.", sep="\t", file=sys.stderr)	
#	print("", sep="\t", file=sys.stderr)
	print("Pairwise SNP Distances:", sep="\t", file=sys.stderr)

	headers = ["Genome Position", "FILTER1", "FILTER", "REF", "ALT", "Gene", "Rv", "Variant Type", "Nucl", "AA", "Score", 'GT', args.vcf1, 'GT', args.vcf2]

	vcf1_sites = fetch_variants(args.vcf1, args.site)
	vcf2_sites = fetch_variants(args.vcf2, args.site)
	
	vcf1_count = 0
	vcf1_only = 0
	for site in vcf1_sites.keys():
		if site not in vcf2_sites.keys():
#			print(vcf1_sites)
			vcf1_only += 1
		vcf1_count += 1

	vcf2_count = 0
	vcf2_only = 0
	for site in vcf2_sites.keys():
		if site not in vcf1_sites.keys():
#			print(vcf2_sites[site])
			vcf2_only += 1
		vcf2_count += 1

	print("VCF1 total", vcf1_count, file=sys.stderr)
	print("VCF1 only", vcf1_only, file=sys.stderr)
	print("VCF2 total", vcf2_count, file=sys.stderr)
	print("VCF2 only", vcf2_only, file=sys.stderr)

	print(*headers, sep="\t")
	
	# define output list positions
	pPos     = 0
	pFilters = 1
	pRef     = 2
	pAlt     = 3
	pGene    = 4
	pRv      = 5
	pVariant = 6
	pNucl    = 7
	pAA      = 8
	pScore   = 9
	pGT      = 10
	pData    = 11

#		output = [record.POS, filters, record.REF, record.ALT[0], gene, rv_name, variant_type, variant_nucl, variant_prot, score]
#		headers = ["Genome Position", "FILTER1", "FILTER", "REF", "ALT", "Gene", "Rv", "Variant Type", "Nucl", "AA", "Score", 'GT', args.vcf1, 'GT', args.vcf2]


	for i in range(0,5000000):
		included = 0
		output = [None] * 13

		if i in vcf1_sites.keys():
			included = 1
#			print("T1")
			
#			output[0:2] = vcf1_sites[i][0:2]
#			output[3:12] = vcf1_sites[i][2:12]

			output[0:2]  = [vcf1_sites[i][j] for j in [pPos, pFilters] ]
			output[3:12] = [vcf1_sites[i][j] for j in [pRef, pAlt, pGene, pRv, pVariant, pNucl, pAA, pScore, pGT, pData] ]

			if i in vcf2_sites.keys():
#				print("T2")
#				output[2] = vcf2_sites[i][1]
#				output[13:14] = vcf2_sites[i][10:12]

				output[2]     = vcf2_sites[i][pFilters] 
				output[13:14] = [vcf2_sites[i][j] for j in [pGT, pData] ]

				if output[5] == '':
					output[4] = vcf2_sites[i][pAlt]
					output[5] = vcf2_sites[i][pGene]
					output[6] = vcf2_sites[i][pRv]
					output[7] = vcf2_sites[i][pVariant]
					output[8] = vcf2_sites[i][pNucl]
					output[9] = vcf2_sites[i][pAA]
					output[10] = vcf2_sites[i][pScore]

		if included == 0:

			if i in vcf2_sites.keys():
#				print("T3 site", i)
				included = 1
#				output[0] = vcf2_sites[i][0]
#				output[2] = vcf2_sites[i][1]
#				output[3:10] = vcf2_sites[i][2:10]
#				output[13:14] = vcf2_sites[i][10:12]

				output[0]     = vcf2_sites[i][pPos]
				output[2]     = vcf2_sites[i][pFilters]
				output[3:11]  = [vcf2_sites[i][j] for j in [pRef, pAlt, pGene, pRv, pVariant, pNucl, pAA, pScore] ]
				output[13:14] = [vcf2_sites[i][j] for j in [pGT, pData] ]

				if i in vcf1_sites.keys():
#					print("T4")
#					output[1] = vcf1_sites[i][1]

					output[1] = vcf1_sites[i][pFilters]

#					output[12:13] = vcf1_sites[i][10:12]
#				output.append(vcf2_sites[i][10])
#				output.append(vcf2_sites[i][11])
		
		if included:	
			conv = lambda i : i or '' 
			result = [conv(i) for i in output] 
			
			if output[11] != output[13]:
				print(*result, sep="\t")
	

def fetch_variants(vcffile, site):
	print("Processing :", vcffile, file=sys.stderr)

	vcf_reader = vcf.Reader(filename=vcffile)
	samples = vcf_reader.samples

	sites = {}

	for record in vcf_reader:
		
		if site:
			if record.POS < site:
				continue
			elif record.POS == site:
				print("Found site:", site, file=sys.stderr)
			else:
				break

		# process annotation field

		annotations = []
		if 'ANN' in record.INFO.keys():
			annotations = record.INFO['ANN']

		gene = ''
		rv_name = ''
		variant_type = ''
		variant_nucl = ''
		variant_prot = ''
		score = ''

		include = 0

#		if len(annotations) > 1:
#			print("Double annotation", record.POS, annotations)
#			exit(0)

		for annotation in annotations:
			fields = annotation.split("|")
			
			gene         = gene + fields[3]
			rv_name      = rv_name + fields[4]
			variant_type = variant_type + fields[1]
			variant_nucl = variant_nucl + fields[9]
			variant_prot = variant_prot + fields[10]
			score        = score + fields[2]
			
#			if variant_type == 'missense_variant':
#				include = 1


		# process samples
		include = 1
		filters = []
		if len(record.FILTER) > 0:
			filters.append(record.FILTER[0])

		# filter out LowQual sites
		if len(record.FILTER) > 0 and record.FILTER[0] == 'LowQual':
			include = 0

#		print("LENGTH:",  filters, len(filters))
		if len(filters) == 0:
			filters = ''
		
		output = [record.POS, filters, record.REF, record.ALT[0], gene, rv_name, variant_type, variant_nucl, variant_prot, score]
			
#		if args.sample1 is not None and args.sample2 is not None:
		call = record.genotype(samples[0])
			
		output = output + fetch(record, call, record.INFO['DP4'])

		sites[record.POS] = output

	return sites

def junk():
	# filter out sites that are the same in both samples
	if call1['GT'] == call2['GT']:
		include = 0
	if call1['AD'][0] > 0 and call1['AD'][1] > 0:
		filters.append('AD1')
	if call2['AD'][0] > 0 and call2['AD'][1] > 0:
		filters.append('AD2')
			
#			if call1['GT'] == '0/1' or call2['GT'] == '0/1':
#				include = 0
	else:
	
		for sample in samples.keys():
			call = record.genotype(samples[sample])
			output = output + fetch(record, call)
		
	output[1] = "|".join(filters)

	if include == 1:
		print(*output, sep="\t")

def fetch(record, call, dp4):
	output_record = [call['GT'], 'PL=' + str(call['PL']) + ', DP=' + str(call['DP']) + ', DP4=' + str(dp4) + ', ADF=' + str(call['ADF']) + ', ADR=' + str(call['ADR']) + ', AD=' + str(call['AD'])]
	return output_record
			

if __name__ == '__main__':
	main()



