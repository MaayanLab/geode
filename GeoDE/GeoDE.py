import warnings
import numpy as np
from sklearn.decomposition import PCA

warnings.filterwarnings("ignore", category=DeprecationWarning) 
warnings.filterwarnings("ignore", category=RuntimeWarning) 


def chdir(data, sampleclass, genes, gamma=1., sort=True):
	"""
	Calculate the characteristic direction for a gene expression dataset
	
	Input:
		data: numpy.array, is the data matrix of gene expression where rows correspond to genes and columns correspond to samples
		sampleclass: list or numpy.array, labels of the samples, it has to be consist of 1 and 2, with 1 being control and 2 being perturbation
				example: sampleclass = [1,1,1,2,2,2]
		genes: list or numpy.array, row labels for genes 
		gamma: float, regulaized term. A parameter that smooths the covariance matrix and reduces potential noise in the dataset
		sort: bool, whether to sort the output by the absolute value of chdir
	Output:
		A list of tuples sorted by the absolute value in descending order characteristic directions of genes.	
	"""
	
	## check input
	data.astype(float)
	sampleclass = np.array(map(int, sampleclass))
	# masks
	m1 = sampleclass == 1 
	m2 = sampleclass == 2

	if type(gamma) not in [float, int]:
		raise ValueError("gamma has to be a numeric number")
	if set(sampleclass) != set([1,2]):
		raise ValueError("sampleclass has to be a list whose elements are in only 1 or 2")
	if m1.sum()<2 or m2.sum()<2:
		raise ValueError("Too few samples to calculate characteristic directions")
	if len(genes) != data.shape[0]:
		raise ValueError("Number of genes does not match the demension of the expression matrix")

	## start to compute
	n1 = m1.sum() # number of controls
	n2 = m2.sum() # number of experiments

	## the difference between experiment mean vector and control mean vector.
	meanvec = data[:,m2].mean(axis=1) - data[:,m1].mean(axis=1) 

	## initialize the pca object
	pca = PCA(n_components=None)
	pca.fit(data.T)

	## compute the number of PCs to keep
	cumsum = pca.explained_variance_ratio_ # explained variance of each PC
	keepPC = len(cumsum[cumsum > 0.001]) # number of PCs to keep

	v = pca.components_[0:keepPC].T # rotated data 
	r = pca.transform(data.T)[:,0:keepPC] # transformed data

	dd = ( np.dot(r[m1].T,r[m1]) + np.dot(r[m2].T,r[m2]) ) / float(n1+n2-2)
	sigma = np.mean(np.diag(dd)) # the scalar covariance

	shrunkMats = np.linalg.inv(gamma*dd + sigma*(1-gamma)*np.eye(keepPC))

	b = np.dot(np.dot(np.dot(v,shrunkMats), v.T), meanvec)
	b /= np.linalg.norm(b) # normalize b to unit vector

	grouped = zip([abs(item) for item in b],b,genes)
	if sort:
		grouped = sorted(grouped,key=lambda x: x[0], reverse=True)

	# return sorted b and genes.
	res = [(item[1],item[2]) for item in grouped]
	return res



def paea(chdir, gmtline, case_sensitive=False):
	"""
	Perform principal angle enrichment analysis (PAEA)
	Input:
		chdir, list of tuples: A characteristic direction returned from chdir function
		gmtline: A list of genes from a gene set 
	Output:
		a list of tuples in the format of:
			the principal angle, the p value
	"""
	if not case_sensitive:
		genes_measured = [gene.upper() for b, gene in chdir]
		gmtline = [gene.upper() for gene in gmtline]
	else: # case sensitive
		genes_measured = [gene for b, gene in chdir]

	if len(set(genes_measured)) != len(chdir):
		raise ValueError('There are duplicated genes in the input genes')
	
	mask = np.in1d(genes_measured, gmtline) # gpos in the R script
	mm = np.where(mask==True)[0]
	m = mask.sum() # number of overlaping genes
	n = len(genes_measured)

	if m > 1 and m < n: # if there is overlap between gene set and genes in chdir
		gsa = np.zeros((n, m)) # Qc in the paper
		for i in range(m):
			gsa[mm[i], i] = 1.

		qb = np.array([b for b,gene in chdir])
		qbqc = np.matrix( np.dot(qb, gsa) ) # Qb.T Qc in paper
		principal_angle = np.linalg.svd(qbqc,compute_uv=False)[0]
		theta = np.arccos(principal_angle)

		# calculation of the null-distribution of principal angles
		m = float(m)
		n = float(n)
		pac = lambda theta: 2.*(1./np.sqrt(2*np.pi))*np.exp((n/2.)*np.log(n/(n - m))+(m/2.)*np.log((n - m)/m)+(1/2.)*np.log(m/(2.*n)*(n - m))+(n-m-1)*np.log(np.sin(theta))+(m-1)*np.log(np.cos(theta)))
		integration_range = np.linspace(0, theta, num=10000, endpoint=True) ## num seems to matter a lot
		p_val = np.trapz(pac(integration_range), integration_range)
		
	else:
		principal_angle = 0.
		p_val = 1.

	return principal_angle, p_val


def paea_wrapper(chdir, gmt_fn, case_sensitive=False):
	"""
	A wrapper function for PAEA gene-set enrichment analysis

	Input:
		chdir: characteristic directions computed by chdir function
		gmt_fn: file name of a gene-set library in GMT format
		case_sensitive: whether gene symbols should be considered as case_sensitive
	Output:
		a sorted list of tuples (term, p_val)
	"""
	## check input:
	if not gmt_fn.endswith('.gmt'):
		raise IOError("The gene-set library file is not in GMT format")

	## read gmt into a dict:
	res = []
	with open (gmt_fn) as f:
		for line in f:
			sl = line.strip().split('\t')
			term = sl[0]
			genes = sl[2:]
			principal_angle, p_val = paea(chdir, genes, case_sensitive=case_sensitive)
			res.append( (term, p_val) )
	## sort terms based on p values in ascending order
	res = sorted(res, key=lambda x:x[1])
	return res

