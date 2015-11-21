# Input files
READ1 = READ1.txt.gz
READ2 = READ2.txt.gz
GENOME = genome.fasta
GBK = genome.gbk
TARGET = target.fasta

# Directories and parameters
FASTQC = $(SOFTDIR)FastQC/fastqc
PICARD = $(SOFTDIR)picard-tools-1.119
GATK = $(SOFTDIR)GATK
PARSNP = $(SOFTDIR)Parsnp-Linux64-v1.2
JAVA7 = $(SOFTDIR)jre1.7.0_76/bin/java
JAVAMEM = 24
CPU = 10
TARGETDEPTH = 10
PLOIDY = 1
THETA = 0.05
SPECIES = ecoli
MAXCOVERAGE = 100
SEED = 100
FILTER = -f "DP > $(TARGETDEPTH)" -g "GQ > 20"

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

SUBSAMPLED1 = $(addsuffix .sub.fq.gz, $(basename $(notdir $(READ1)) .txt.gz))
SUBSAMPLED2 = $(addsuffix .sub.fq.gz, $(basename $(notdir $(READ2)) .txt.gz))

$(SUBSAMPLED1): $(READ1) $(READ2) $(GENOME)
	sample=$$($(SRCDIR)/get_subsample $(GENOME) $$(interleave_pairs $(READ1) $(READ2) | count_seqs | awk '{print $$2}') --coverage 100) && \
        seqtk sample -s$(SEED) $(READ1) $$sample > $(SUBSAMPLED1) && \
	seqtk sample -s$(SEED) $(READ2) $$sample > $(SUBSAMPLED2)	

# Alignment
GINDEX = $(GENOME).bwt
$(GINDEX): $(GENOME)
	bwa index $(GENOME)

ALIGNMENT = aln.sam
$(ALIGNMENT): $(GINDEX) $(SUBSAMPLED1)
	bwa mem -t $(CPU) $(GENOME) $(SUBSAMPLED1) $(SUBSAMPLED2) > $(ALIGNMENT)

SORTEDALIGN = aln.sorted.bam
$(SORTEDALIGN): $(ALIGNMENT)
	samtools view -bS $(ALIGNMENT) -q 25 -f 2 -F 256 -o aln.bam && \
	samtools sort aln.bam aln.sorted

BINDEX = $(SORTEDALIGN).bai
$(BINDEX): $(SORTEDALIGN)
	samtools index $(SORTEDALIGN)

DEDUPALIGN = aln.dedup.bam
$(DEDUPALIGN): $(SORTEDALIGN) $(BINDEX)
	java -Xmx4g -jar $(PICARD)/MarkDuplicates.jar INPUT=$(SORTEDALIGN) OUTPUT=$(DEDUPALIGN) METRICS_FILE=metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT

DINDEX = $(DEDUPALIGN).bai
$(DINDEX): $(DEDUPALIGN)
	samtools index $(DEDUPALIGN)

REALIGN = realn.dedup.bam
$(REALIGN): $(DEDUPALIGN) $(GENOME) $(DINDEX)
	java -Xmx$(JAVAMEM)g -jar $(PICARD)/AddOrReplaceReadGroups.jar \
	   I=$(DEDUPALIGN) O=$(DEDUPALIGN).group.bam RGPL=illumina RGLB=foo RGPU=run \
	   RGSM=anysample CREATE_INDEX=true && \
	java -jar $(PICARD)/CreateSequenceDictionary.jar R=$(GENOME) \
	    O=$(basename $(GENOME)).dict GENOME_ASSEMBLY=genome SPECIES=$(SPECIES)
	samtools faidx $(GENOME)
	$(JAVA7) -Xmx$(JAVAMEM)g -jar $(GATK)/GenomeAnalysisTK.jar \
	   -T RealignerTargetCreator \
	   -R $(GENOME) \
	   -I $(DEDUPALIGN).group.bam \
	   -o $(DEDUPALIGN).group.bam.intervals && \
	$(JAVA7) -Xmx$(JAVAMEM)g -jar $(GATK)/GenomeAnalysisTK.jar \
	   -T IndelRealigner \
	   -R $(GENOME) \
	   -I $(DEDUPALIGN).group.bam \
	   -targetIntervals $(DEDUPALIGN).group.bam.intervals \
	   -o $(DEDUPALIGN).realn.bam && \
	samtools calmd -brA $(DEDUPALIGN).realn.bam $(GENOME) > $(REALIGN)

RINDEX = $(REALIGN).bai
$(RINDEX): $(REALIGN)
	samtools index $(REALIGN)

MAPVARIANTS = map.vcf
$(MAPVARIANTS): $(RINDEX) $(REALIGN) $(GENOME)
	freebayes -f $(GENOME) --ploidy $(PLOIDY) --theta $(THETA) --genotype-qualities --standard-filters $(REALIGN) > raw.vcf
	vcffilter $(FILTER) raw.vcf > $(MAPVARIANTS) 
map: $(MAPVARIANTS)

TINDEX = $(TARGET).bwt
$(TINDEX): $(TARGET)
	bwa index $(TARGET)

TALIGNMENT = aln.target.sam
$(TALIGNMENT): $(TINDEX) $(SUBSAMPLED1)
	bwa mem -t $(CPU) $(TARGET) $(SUBSAMPLED1) $(SUBSAMPLED2) > $(TALIGNMENT)

TSORTEDALIGN = aln.target.sorted.bam
$(TSORTEDALIGN): $(TALIGNMENT)
	samtools view -bS $(TALIGNMENT) -q 25 -f 2 -F 256 -o aln.target.bam && \
	samtools sort aln.target.bam aln.target.sorted

BTINDEX = $(TSORTEDALIGN).bai
$(BTINDEX): $(TSORTEDALIGN)
	samtools index $(TSORTEDALIGN)

MASK = mask.bed
$(MASK): $(TSORTEDALIGN)
	bedtools genomecov -ibam $(TSORTEDALIGN) -bg | awk '{if ($$4 < $(TARGETDEPTH)) print $$0}' > $(MASK)

REPEATS = repeats.bed
$(REPEATS): $(GENOME)
	nucmer --maxmatch --nosimplify $(GENOME) $(GENOME) && \
	show-coords -r -T out.delta -H | awk '{if ($$1 != $$3 && $$2 != $$4) print $$0}' > repeats.txt && \
	awk '{print $$8"\t"$$1"\t"$$2}' repeats.txt > $(REPEATS)

PARSNPOUT = parsnp/parsnp.ggr
$(PARSNPOUT): $(GENOME) $(TARGET) $(MASK) $(REPEATS)
	mkdir -p genomes && \
	bedtools maskfasta -fi $(GENOME) -bed $(REPEATS) -fo genomes/$(shell basename $(GENOME))
	bedtools maskfasta -fi $(TARGET) -bed $(MASK) -fo genomes/$(shell basename $(TARGET))
	$(PARSNP)/parsnp -r genomes/$(GENOME) -d genomes -p $(CPU) -v -c -o parsnp

ALIGNVARIANTS = align.vcf
$(ALIGNVARIANTS): $(PARSNPOUT) 
	harvesttools -i $(PARSNPOUT) -V $(ALIGNVARIANTS).vcf
	$(SRCDIR)/parsnp2vcf $(ALIGNVARIANTS).vcf $(ALIGNVARIANTS) --template $(ALIGNVARIANTS).vcf
align: $(ALIGNVARIANTS)

PARSNPOUT1 = parsnp1/parsnp.ggr
$(PARSNPOUT1): $(GENOME) $(TARGET) $(REPEATS)
	mkdir -p genomes && \
	bedtools maskfasta -fi $(GENOME) -bed $(REPEATS) -fo genomes/$(shell basename $(GENOME))
	cp $(TARGET) genomes
	$(PARSNP)/parsnp -r genomes/$(GENOME) -d genomes -p $(CPU) -v -c -o parsnp1

ALIGNVARIANTS1 = align.nomap.vcf
$(ALIGNVARIANTS1): $(PARSNPOUT1) 
	harvesttools -i $(PARSNPOUT1) -V $(ALIGNVARIANTS1).vcf
	$(SRCDIR)/parsnp2vcf $(ALIGNVARIANTS1).vcf $(ALIGNVARIANTS1) --template $(ALIGNVARIANTS1).vcf
alignnoreads: $(ALIGNVARIANTS1)

BRESEQOUT = $(CURDIR)/output/output.gd
BRESEQVARIANTS = breseq.vcf
PREAD1 = PREAD1.txt.gz
$(PREAD1): $(READ1)
	zcat $(READ1) > $(PREAD1)
PREAD2 = PREAD2.txt.gz
$(PREAD2): $(READ2)
	zcat $(READ2) > $(PREAD2)

$(BRESEQOUT): $(PREAD1) $(PREAD2) $(GBK)
	breseq -r $(GBK) $(PREAD1) $(PREAD2) -j $(CPU)

$(BRESEQVARIANTS): $(BRESEQOUT) $(GBK)
	gdtools GD2VCF -r $(GBK) $(BRESEQOUT) -o $(BRESEQVARIANTS)
breseq: $(BRESEQVARIANTS)

all: fastqc trim map align alignnoreads breseq

.PHONY: all fastqc trim map align alignnoreads breseq
