Scripts by Robert Williamson (http://www.genomicconflict.com/wiki/index.php?title=Category:Programs and https://github.com/fabbyrob/science/tree/master/pileup_analyzers):

annotation.py
bootstrapRegions.py
filterSummary.py	only works in Python v2.7 and together with summary.py
summary.py
vcf.py	vcf parser by R. Williamson
vcfSummarizer.py	modified to take gzipped files, and to filter for 'RGQ', the genotype quality of non-variant sites.


directory PyVCF --> PyVCF scripts (vcf parser, http://pyvcf.readthedocs.org):

cparse.pyx
filters.py
__init__.py
model.py
sample_filter.py
utils.py
vcf.py	original name on PyVCF github website: parser.py. Renamed to be recognized by R. Williamsons scripts.
