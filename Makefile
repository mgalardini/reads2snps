# Input files
READ1 = READ1.txt.gz
READ2 = READ2.txt.gz
GENOME = genome.fasta

# Directories and parameters
FASTQC = FastQC/fastqc
PLOIDY = 1
THETA = 0.05
MAXCOVERAGE = 100
SEED = 100

# Anything below this point should not be changed

# Directories
SRCDIR = $(CURDIR)/src

QCDIR = $(CURDIR)/QC
$(QCDIR):
	mkdir -p $(QCDIR)

TRIMDIR = $(CURDIR)/trimmed
$(TRIMDIR):
	mkdir -p $(TRIMDIR)

# QC
QCREAD1 = $(QCDIR)/$(addsuffix _fastqc.zip, $(basename $(notdir $(READ1))))
QCREAD2 = $(QCDIR)/$(addsuffix _fastqc.zip, $(basename $(notdir $(READ2))))

$(QCREAD1): $(QCDIR) $(READ1)
	$(FASTQC) --outdir $(QCDIR) $(READ1)
$(QCREAD2): $(QCDIR) $(READ2)
	$(FASTQC) --outdir $(QCDIR) $(READ2)
fastqc: $(QCREAD1) $(QCREAD2)

# Trimming
TREAD1 = $(addsuffix .fq.gz, $(TRIMDIR)/$(basename $(notdir $(READ1)) .txt))
TREAD2 = $(addsuffix .fq.gz, $(TRIMDIR)/$(basename $(notdir $(READ2)) .txt))

$(TREAD1): $(TRIMDIR) $(READ1) $(READ2)
	interleave_pairs $(READ1) $(READ2) | \
	trim_edges -l 9 --paired_reads | \
	trim_quality -q 20 -w 5 --paired_reads | \
	deinterleave_pairs -z -o $(TREAD1) $(TREAD2)
trim: $(TREAD1)

SUBSAMPLED1 = $(addsuffix .sub.fq.gz, $(TRIMDIR)/$(basename $(notdir $(READ1)) .txt))
SUBSAMPLED2 = $(addsuffix .sub.fq.gz, $(TRIMDIR)/$(basename $(notdir $(READ2)) .txt))

$(SUBSAMPLED1): $(TREAD1) $(GENOME)
	sample=$$($(SRCDIR)/get_subsample $(GENOME) $$(interleave_pairs $(TREAD1) $(TREAD2) | count_seqs | awk '{print $$2}') --coverage 100) && \
        seqtk sample -s$(SEED) $(TREAD1) $$sample > $(SUBSAMPLED1) && \
	seqtk sample -s$(SEED) $(TREAD2) $$sample > $(SUBSAMPLED2)	

# Alignment
GINDEX = $(GENOME).bwt
$(GINDEX): $(GENOME)
	bwa index $(GENOME)

ALIGNMENT = aln.sam
$(ALIGNMENT): $(GINDEX) $(SUBSAMPLED1)
	bwa mem $(GENOME) $(SUBSAMPLED1) $(SUBSAMPLED2) > $(ALIGNMENT)

SORTEDALIGN = aln.sorted.bam
$(SORTEDALIGN): $(ALIGNMENT)
	samtools view -bS $(ALIGNMENT) -o aln.bam && \
	samtools sort aln.bam aln.sorted

BINDEX = $(SORTEDALIGN).bai
$(BINDEX): $(SORTEDALIGN)
	samtools index $(SORTEDALIGN)

VARIANTS = var.vcf
$(VARIANTS): $(BINDEX) $(SORTEDALIGN) $(GENOME)
	freebayes -f $(GENOME) --ploidy $(PLOIDY) --theta $(THETA) $(SORTEDALIGN) > $(VARIANTS)
variants: $(VARIANTS)

all: fastqc trim variants

.PHONY: all fastqc trim variants
