from geode import *
from pprint import pprint
import numpy as np
import os.path

## read expression data
mat = []
genes = []
with open (os.path.join(os.path.dirname(__file__), 'example_expression_data.txt')) as f:
	next(f)
	for line in f:
		sl = line.strip().split('\t')
		gene = sl[0]
		row = list(map(float, sl[1:]))
		genes.append(gene)
		mat.append(row)
mat = np.array(mat)
col_labels = ['1','1','1','2','2','2']

## compute characteristic direction
chdir_res = chdir(mat, col_labels, genes, calculate_sig=0, nnull=100)
pprint(chdir_res)
# perform PAEA gene-set enrichment analysis
paea_res = paea_wrapper(chdir_res, 'GeneOntology_BP.gmt')
# look at the top enriched GO terms
pprint(paea_res[0:10])
