# ChIP-seq normalisation using spike-in

ifdef SLURM_JOB_ID
SHELL := srun
.SHELLFLAGS := -N 1 -n 1 --cpus-per-task=4  bash -c
else
SHELL := /bin/bash
endif

nthread := 10

.PHONY: all test clean

genome_path = $(HOME)/data/genome

target_genome := ce10
spikein_genome := cb3
concatenated_genome = $(target_genome)-$(spikein_genome)

target_genome_fa = $(call GET_GENOME_FA,$(target_genome))
spikein_genome_fa = $(call GET_GENOME_FA,$(spikein_genome))
concatenated_genome_fa = $(call GET_GENOME_FA,$(concatenated_genome))

target_genome_size = $(call GET_GENOME_SIZE,$(target_genome))
spikein_genome_size = $(call GET_GENOME_SIZE,$(spikein_genome))
concatenated_genome_size = $(call GET_GENOME_SIZE,$(concatenated_genome))

blacklists = $(foreach g,$(target_genome),$(spikein_genome),$(g).blacklist.bed)

sample_tbl := samples.txt
samples = $(shell bioawk -t -c hdr 'NR>1{print $$sample}' $(sample_tbl))
input_samples = $(shell bioawk -t -c hdr '$$factor=="input"{print $$sample}' $(sample_tbl))
chip_samples = $(shell bioawk -t -c hdr 'NR>1 && $$factor!="input"{print $$sample}' $(sample_tbl))
factors = $(shell bioawk -t -c hdr 'NR>1{print $$factor}' $(sample_tbl) | sort | uniq | grep -vi input)

aligned_bams = $(foreach s,$(samples),$(s)/$(s).aligned.bam)
target_bams = $(foreach s,$(samples),$(s)/$(s).aligned.$(target_genome).bam)
spikein_bams = $(foreach s,$(samples),$(s)/$(s).aligned.$(spikein_genome).bam)
split_bams = $(target_bams) $(spikein_bams)

target_input_bw = $(call GET_POOLED_INPUT,$(target_genome))
spikein_input_bw = $(call GET_POOLED_INPUT,$(spikein_genome))

peaks = $(foreach f,$(factors),$(foreach g,$(target_genome) $(spikein_genome),$(f).$(g).peaks.bed))
spikein_peak_heights = $(foreach s,$(chip_samples),$(s)/$(s).$(spikein_genome).mq$(MQ).avgPeakHeight)

total_counts = $(foreach s,$(samples),$(foreach g,$(target_genome) $(spikein_genome),$(s)/$(s).$(g).mq$(MQ).total.readCount))
blacklist_counts = $(subst total,blacklist,$(total_counts))
peak_counts = $(foreach s,$(chip_samples),$(foreach g,$(target_genome) $(spikein_genome),$(s)/$(s).$(g).mq$(MQ).peak.readCount))

stats_tbl = norm_stats.txt

norm_tracks = $(foreach s,$(chip_samples),$(foreach n,$(NORM_METHOD),$(s)/$(s).mq$(MQ).spknorm$(n).bw))

# macros
GET_GENOME_FA = $(genome_path)/$(1)/$(1).fa
GET_GENOME_SIZE = $(genome_path)/$(1)/$(1).fa.fai
GET_POOLED_INPUT = $(1).pooled_input.bw
GET_SAMPLE_FOR_FACTOR = $(shell bioawk -t -c hdr '$$factor=="$(1)" {print $$sample}' $(sample_tbl))
GET_FACTOR_FOR_SAMPLE = $(shell bioawk -t -c hdr '$$sample=="$(1)" {print $$factor}' $(sample_tbl))
GET_INPUT_FOR_CHIP = $(shell bioawk -t -c hdr '$$sample=="$(1)" {print $$input}' $(sample_tbl))

# parameters
MQ := 10
EXTSIZE := 200
MACS2_PEAK_C := 2
NORM_METHOD := 0 1 2 3 4 5 6 7

all: $(stats_tbl) $(norm_tracks)

clean:
	-rm -f $(stats_tbl)

test:
	@echo $(split_bams)

.PRECIOUS: $(split_bams)

.INTERMEDIATE: $(foreach s,$(samples),$(s)/$(s).sortp.bam) $(foreach s,$(samples),$(s)/$(s).bam) $(foreach s,$(samples),$(s)/$(s).sai)

define ALIGN
$(1)/$(1).sai: $$(concatenated_genome_fa) fastq/$(1).fastq.gz
	bwa aln -t $$(nthread) $$^ > $$@

$(1)/$(1).bam: $$(concatenated_genome_fa) $(1)/$(1).sai fastq/$(1).fastq.gz
	bwa samse $$^ | samtools view -b -o $$@ - && rm -f $(1)/$(1).sai

$(1)/$(1).sortp.bam: $(1)/$(1).bam
	sambamba sort -o $$@ $$< && rm -f $$<

$(1)/$(1).aligned.bam: $(1)/$(1).sortp.bam
	sambamba markdup --hash-table-size=1000000 --overflow-list-size=1000000 $$< $$@ && rm -f $$<
endef

$(foreach s,$(samples),$(eval $(call ALIGN,$(s))))

define SPLIT
$(1)/$(1).aligned.$(2).bam: $(1)/$(1).aligned.bam
	grep -P '^$(2)' $$(concatenated_genome_size) | cut -f1 | parallel -k --xargs sambamba view -t 2 $$< {} \
		| sed 's/$(2)_chr/chr/g' | samtools view -@ 4 -b -t $$(call GET_GENOME_SIZE,$(2)) -o $$@ - \
		&& sambamba index $$@
endef

$(foreach s,$(samples),$(foreach g,$(target_genome) $(spikein_genome),$(eval $(call SPLIT,$(s),$(g)))))

define PILEUP
$(1)/$(1).aligned.$(2).macs2_ext$$(EXTSIZE).bdg: $(1)/$(1).aligned.$(2).bam
	macs2 pileup -i $$< -o $$@ --extsize $$(EXTSIZE)

$(1)/$(1).aligned.$(2).macs2_ext$$(EXTSIZE).bw: $(1)/$(1).aligned.$(2).macs2_ext$$(EXTSIZE).bdg
	bedGraphToBigWig $$< $$(call GET_GENOME_SIZE,$(2)) $$@
endef

$(foreach s,$(samples),$(foreach g,$(target_genome) $(spikein_genome),$(eval $(call PILEUP,$(s),$(g)))))

define DIVIDE_BY_INPUT
$(1)/$(1).aligned.$(2).macs2_ext$$(EXTSIZE).input_norm.bw: $(1)/$(1).aligned.$(2).macs2_ext$$(EXTSIZE).bw $$(call GET_POOLED_INPUT,$(2))
	bwDivide -m 0.25 -o $$@ $$^

$(1)/$(1).aligned.$(2).macs2_ext$$(EXTSIZE).input_norm.median_norm.bw: $(1)/$(1).aligned.$(2).macs2_ext$$(EXTSIZE).bw $$(call GET_POOLED_INPUT,$(2))
	bwDivide -M -m 0.25 -o $$@ $$^
endef

$(foreach s,$(samples),$(foreach g,$(target_genome) $(spikein_genome),$(eval $(call DIVIDE_BY_INPUT,$(s),$(g)))))

define CALL_PEAK
$(1)/$(1).aligned.$(2).macs2_ext$$(EXTSIZE).input_norm.median_norm.macs2_c$$(MACS2_PEAK_C).narrowPeak: $(1)/$(1).aligned.$(2).macs2_ext$$(EXTSIZE).input_norm.median_norm.bw
	bigWigToBedGraph $$< $$(subst .bw,.bdg,$$<)
	macs2 bdgpeakcall -i $$(subst .bw,.bdg,$$<) -c $$(MACS2_PEAK_C) -o $$@ && rm -f $$(subst .bw,.bdg,$$<)
endef

$(foreach s,$(chip_samples),$(foreach g,$(target_genome) $(spikein_genome),$(eval $(call CALL_PEAK,$(s),$(g)))))

define MERGE_PEAK
$(1).$(2).peaks.bed: $$(foreach s,$$(call GET_SAMPLE_FOR_FACTOR,$(1)),$$(s)/$$(s).aligned.$(2).macs2_ext$$(EXTSIZE).input_norm.median_norm.macs2_c$$(MACS2_PEAK_C).narrowPeak)
	cat $$^ | sort -k1,1 -k2,2n -k3,3n | bedtools merge | sort -k1,1 -k2,2n -k3,3n | bed3to6 | cut -f1-4 > $$@
endef

$(foreach f,$(factors),$(foreach g,$(target_genome) $(spikein_genome),$(eval $(call MERGE_PEAK,$(f),$(g)))))

define COUNT_ALL
$(1)/$(1).$(2).mq$$(MQ).total.readCount: $(1)/$(1).aligned.$(2).bam
	samtools view -c -q $$(MQ) $$< > $$@
endef

$(foreach s,$(samples),$(foreach g,$(target_genome) $(spikein_genome),$(eval $(call COUNT_ALL,$(s),$(g)))))

define COUNT_BLACKLIST
$(1)/$(1).$(2).mq$$(MQ).blacklist.readCount: $(1)/$(1).aligned.$(2).bam $(2).blacklist.bed
	cut -f1-3 $(2).blacklist.bed | awk '{print $$$$1":"$$$$2"-"$$$$3}' | parallel -k --xargs samtools view -@ 2 -c -q $$(MQ) $$< {} > $$@
endef

$(foreach s,$(samples),$(foreach g,$(target_genome) $(spikein_genome),$(eval $(call COUNT_BLACKLIST,$(s),$(g)))))

define COUNT_CHIP_PEAK
$(1)/$(1).$(2).mq$$(MQ).peak.readCount: $(1)/$(1).aligned.$(2).bam $$(call GET_FACTOR_FOR_SAMPLE,$(1)).$(2).peaks.bed
	cut -f1-3 $$(call GET_FACTOR_FOR_SAMPLE,$(1)).$(2).peaks.bed | awk '{print $$$$1":"$$$$2"-"$$$$3}' \
		| parallel -k --xargs samtools view -@ 2 -c -q $$(MQ) $$< {} > $$@
endef

$(foreach s,$(chip_samples),$(foreach g,$(target_genome) $(spikein_genome),$(eval $(call COUNT_CHIP_PEAK,$(s),$(g)))))

define CALC_AVG_PEAK_HEIGHT
$(1)/$(1).$(2).mq$$(MQ).avgPeakHeight: $(1)/$(1).aligned.$(2).macs2_ext$$(EXTSIZE).input_norm.median_norm.bw $$(call GET_FACTOR_FOR_SAMPLE,$(1)).$(2).peaks.bed
	bigWigAverageOverBed -sampleAroundCenter=50 $$^ /dev/stdout | bioawk -t '{s=s+$$$$6}END{print s/(NR-1)}' > $$@
endef

$(foreach s,$(chip_samples),$(eval $(call CALC_AVG_PEAK_HEIGHT,$(s),$(spikein_genome))))

$(stats_tbl): $(sample_tbl) $(total_counts) $(blacklist_counts) $(peak_counts) $(spikein_peak_heights)
	python calcNormFactor.py -t $(target_genome) -s $(spikein_genome) -q $(MQ) $< > $@

define GET_NORM_STATS
$(1)/$(1).mq$$(MQ).normStats: $(stats_tbl)
	bioawk -t -c hdr 'NR==1 || $$$$sample=="$(1)"' $$< > $$@
endef

$(foreach s,$(chip_samples),$(eval $(call GET_NORM_STATS,$(s))))

define GET_NORM_TRACK
$(1)/$(1).mq$$(MQ).spknorm0.bw: $(1)/$(1).aligned.$$(target_genome).macs2_ext$$(EXTSIZE).input_norm.bw $(1)/$(1).mq$$(MQ).normStats
	bigWigToBedGraph $$< /dev/stdout \
		| bioawk -t -v f=`bioawk -t -c hdr 'NR>1{print $$$$norm_factor}' $(1)/$(1).mq$$(MQ).normStats` '{print $$$$1,$$$$2,$$$$3,$$$$4*f}' > $$(subst .bw,.bdg,$$@)
	bedGraphToBigWig $$(subst .bw,.bdg,$$@) $$(call GET_GENOME_SIZE,$$(target_genome)) $$@ && rm -f $$(subst .bw,.bdg,$$@)

$(1)/$(1).mq$$(MQ).spknorm1.bw: $(1)/$(1).aligned.$$(target_genome).macs2_ext$$(EXTSIZE).input_norm.bw $(1)/$(1).mq$$(MQ).normStats
	bigWigToBedGraph $$< /dev/stdout \
		| bioawk -t -v f=`bioawk -t -c hdr 'NR>1{print $$$$norm_factor_rmblack}' $(1)/$(1).mq$$(MQ).normStats` '{print $$$$1,$$$$2,$$$$3,$$$$4*f}' > $$(subst .bw,.bdg,$$@)
	bedGraphToBigWig $$(subst .bw,.bdg,$$@) $$(call GET_GENOME_SIZE,$$(target_genome)) $$@ && rm -f $$(subst .bw,.bdg,$$@)

$(1)/$(1).mq$$(MQ).spknorm2.bw: $(1)/$(1).aligned.$$(target_genome).macs2_ext$$(EXTSIZE).input_norm.bw $(1)/$(1).mq$$(MQ).normStats
	bigWigToBedGraph $$< /dev/stdout \
		| bioawk -t -v f=`bioawk -t -c hdr 'NR>1{print $$$$norm_factor_peak}' $(1)/$(1).mq$$(MQ).normStats` '{print $$$$1,$$$$2,$$$$3,$$$$4*f}' > $$(subst .bw,.bdg,$$@)
	bedGraphToBigWig $$(subst .bw,.bdg,$$@) $$(call GET_GENOME_SIZE,$$(target_genome)) $$@ && rm -f $$(subst .bw,.bdg,$$@)

$(1)/$(1).mq$$(MQ).spknorm3.bw: $(1)/$(1).aligned.$$(target_genome).macs2_ext$$(EXTSIZE).input_norm.bw $(1)/$(1).mq$$(MQ).normStats
	bigWigToBedGraph $$< /dev/stdout \
		| bioawk -t -v f=`bioawk -t -c hdr 'NR>1{print $$$$norm_factor_bg}' $(1)/$(1).mq$$(MQ).normStats` '{print $$$$1,$$$$2,$$$$3,$$$$4*f}' > $$(subst .bw,.bdg,$$@)
	bedGraphToBigWig $$(subst .bw,.bdg,$$@) $$(call GET_GENOME_SIZE,$$(target_genome)) $$@ && rm -f $$(subst .bw,.bdg,$$@)

$(1)/$(1).mq$$(MQ).spknorm4.bw: $(1)/$(1).aligned.$$(target_genome).macs2_ext$$(EXTSIZE).input_norm.bw $(1)/$(1).mq$$(MQ).normStats
	bigWigToBedGraph $$< /dev/stdout \
		| bioawk -t -v f=`bioawk -t -c hdr 'NR>1{print $$$$norm_factor_snr}' $(1)/$(1).mq$$(MQ).normStats` '{print $$$$1,$$$$2,$$$$3,$$$$4*f}' > $$(subst .bw,.bdg,$$@)
	bedGraphToBigWig $$(subst .bw,.bdg,$$@) $$(call GET_GENOME_SIZE,$$(target_genome)) $$@ && rm -f $$(subst .bw,.bdg,$$@)

$(1)/$(1).mq$$(MQ).spknorm5.bw: $(1)/$(1).aligned.$$(target_genome).macs2_ext$$(EXTSIZE).input_norm.bw $(1)/$(1).mq$$(MQ).normStats
	bigWigToBedGraph $$< /dev/stdout \
		| bioawk -t -v f=`bioawk -t -c hdr 'NR>1{print $$$$norm_factor_rmblack_snr}' $(1)/$(1).mq$$(MQ).normStats` '{print $$$$1,$$$$2,$$$$3,$$$$4*f}' > $$(subst .bw,.bdg,$$@)
	bedGraphToBigWig $$(subst .bw,.bdg,$$@) $$(call GET_GENOME_SIZE,$$(target_genome)) $$@ && rm -f $$(subst .bw,.bdg,$$@)

$(1)/$(1).mq$$(MQ).spknorm6.bw: $(1)/$(1).aligned.$$(target_genome).macs2_ext$$(EXTSIZE).input_norm.bw $(1)/$(1).mq$$(MQ).normStats
	bigWigToBedGraph $$< /dev/stdout \
		| bioawk -t -v f=`bioawk -t -c hdr 'NR>1{print $$$$norm_factor_peak_snr}' $(1)/$(1).mq$$(MQ).normStats` '{print $$$$1,$$$$2,$$$$3,$$$$4*f}' > $$(subst .bw,.bdg,$$@)
	bedGraphToBigWig $$(subst .bw,.bdg,$$@) $$(call GET_GENOME_SIZE,$$(target_genome)) $$@ && rm -f $$(subst .bw,.bdg,$$@)

$(1)/$(1).mq$$(MQ).spknorm7.bw: $(1)/$(1).aligned.$$(target_genome).macs2_ext$$(EXTSIZE).input_norm.bw $(1)/$(1).mq$$(MQ).normStats
	bigWigToBedGraph $$< /dev/stdout \
		| bioawk -t -v f=`bioawk -t -c hdr 'NR>1{print $$$$norm_factor_bg_snr}' $(1)/$(1).mq$$(MQ).normStats` '{print $$$$1,$$$$2,$$$$3,$$$$4*f}' > $$(subst .bw,.bdg,$$@)
	bedGraphToBigWig $$(subst .bw,.bdg,$$@) $$(call GET_GENOME_SIZE,$$(target_genome)) $$@ && rm -f $$(subst .bw,.bdg,$$@)
endef

$(foreach s,$(chip_samples),$(eval $(call GET_NORM_TRACK,$(s))))
