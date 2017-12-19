Transgene (Translation of Genomic mutations) is a python package to process input mutations in genomic and transcriptomic 
space to produce (2n-1)-mer peptides that can be used to identify neoepitopes of immunotherapeutic potential. Transgene 
can naively handle Single Nucleotide Variants (SNVs) and MultiAllelic Variants (MAVs), and will also handle Short Insertions 
and Deletions (INDELs) and fusion genes if a genomic FASTA and an annotation GTF are provided.

Transgene accepts a vcf of SNVs and INDELs, and a bedpe file containing FUSIONs, and produces output peptide FASTAs 
containing ImmunoActive Regions (IARs), the regions of a protein containing a non-synonymous mutation and the flanking 
WT sequence.  Transgene will also attempt to chain mutations within 3n bases of each other (where n is the length of 
epitope expected to be tested against an MHC binding prediction software). Transgene can currently chain SNVs, including
Multi Nucleotide Polymorphisms (MNPs) that affect a single codon, MAVs, and INDELs in any combination. We plan on chain 
fusions as well in the future.

Transgene uses an optional RNA-Seq BAM file to phase the chained mutations, to ensure that only the "real" combinations of 
mutations are produced. The RNA-Seq file is also used to filter mutants that show low or no rna-seq expression.

Transgene can process a vcf for OxoG variants provided an input Tumor DNA BAM file is provided.

To learn more about how to run transgene, try

         python transgene.py --help
