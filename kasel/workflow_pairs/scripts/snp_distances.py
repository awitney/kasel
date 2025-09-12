
import sys, os
import argparse

from ete3 import Tree

## Collect command line arguments

parser = argparse.ArgumentParser()

parser.add_argument("-t", "--tree", help="Pyjar tree file")
parser.add_argument("-c", "--sitecount", help="Site count in infile matrix", default=1)
parser.add_argument("-s1", "--sample1", help="Sample name 1")
parser.add_argument("-s2", "--sample2", help="Sample name 2")

args = parser.parse_args()

treefile  = args.tree
sitecount = args.sitecount
sample1  = args.sample1
sample2  = args.sample2

## Load tree

t = Tree(treefile, format=1)

## collect leaves

leaves = []

for leaf in t.iter_leaves():
	leaves.append(leaf)

leaves.sort(key=lambda x: x.name)

## Start output

row = ["Name"]

for leaf1 in leaves:
	
	if sample2 is not None:
		if leaf1.name != sample2:
			continue
	
	row.append(leaf1.name)

print("\t".join(row))

## Process leaves

for leaf1 in leaves:

	if sample1 is not None:
		if leaf1.name != sample1:
			continue

	row = []
	row.append(leaf1.name)
	for leaf2 in leaves:
	
		if sample2 is not None:
			if leaf2.name != sample2:
				continue
		
		dist = leaf1.get_distance(leaf2)
		snps = "%.0f" % (float(dist) * int(sitecount))
#		print("%s\t%s\t%.3f\t%.0f" % (leaf1.name, leaf2.name, leaf1.get_distance(leaf2), snps))
#		print("%s\t%s\t%.0f\t%.0f\t%.0f" % (leaf1.name, leaf2.name, leaf1.get_distance(leaf2), float(snps), dist), file=sys.stderr)
		row.append(str(snps))
	
	print("\t".join(row))

