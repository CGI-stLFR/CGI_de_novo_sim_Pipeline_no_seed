# Snakemake Pipeline for sim de novo datasets
# Start from absolute scratch
# So need to make the assembly using the newest software (config file)
# prepend commands with export Sentieon stuff.

configfile: "config.yaml"
shell.prefix('export SENTIEON_INSTALL={}; export SENTIEON_LICENSE={}; '.format(config['params']['sentieon_install'], config['params']['sentieon_license']))


# Rule for run all, could be shortened with a function
# These are the targets we want to make
rule run_all:
    input:
        "Quast_Contigs_0/quast.log",
        "Quast_Contigs_1/quast.log",
        "Contig_Binning/HapA_to_RefB/hapA_refB.vcf",
        "Contig_Binning/HapB_to_RefB/hapB_refB.vcf",
        "Contig_Binning/HapA_to_RefA/hapA_refA.vcf",
        "Contig_Binning/HapB_to_RefA/hapB_refA.vcf",
        "Contig_Binning/Quast_Contigs_HapA/quast.log",
        "Contig_Binning/Quast_Contigs_HapB/quast.log",
        "Scaffolds/HapA_to_RefB/hapA_refB.vcf",
        "Scaffolds/HapB_to_RefB/hapB_refB.vcf",
        "Scaffolds/HapA_to_RefA/hapA_refA.vcf",
        "Scaffolds/HapB_to_RefA/hapB_refA.vcf",
        "Scaffolds/Quast_Scaffolds_HapA/quast.log",
        "Scaffolds/Quast_Scaffolds_HapB/quast.log",
        "Contig_Binning/tmp_contigs_hap0_binA.txt",
        "Contig_Binning/tmp_contigs_hap0_binB.txt",
        "Contig_Binning/tmp_contigs_hap1_binA.txt",
        "Contig_Binning/tmp_contigs_hap1_binB.txt",
        "hap0_phase_n50.txt",
        "hap1_phase_n50.txt",
        "Contig_Binning/SV_Calling/PAF_BED_Alignments/hapA_ref.bed",
        "Contig_Binning/SV_Calling/PAF_BED_Alignments/hapB_ref.bed",
        "Contig_Binning/SV_Calling/BAM_Alignments/hapA_ref.bam",
        "Contig_Binning/SV_Calling/BAM_Alignments/hapB_ref.bam",
        "Contig_Binning/SV_Calling/Stratification_Results/results.png",
        "Contig_Binning/SV_Calling/SV_Comparison_Diploid/summary.txt"


# Rule for read clustering
# mapping reads ot the seed library
rule read_clustering:
    input:
        expand("{sample}", sample=config['samples']['fastq'])
    output:
        "Assembly/cg50X_SentieonIter0.bam"
    log:
        "Assembly/run.log"
    threads:
        config['threads']['assembly']
    params:
        seed_lib = config['params']['seed_lib'],
        readgroup = r'@RG\tID:{0}\tSM:{0}\tPL:{1}'.format(config['samples']['date'],
                                                          config['params']['platform'])
    shell:
        "($SENTIEON_INSTALL/bin/sentieon bwa mem -R '{params.readgroup}' "
            "-M -t {threads} -K 10000000 -Y {params.seed_lib} {input} || echo -n 'BWA error') | "
            "extract | "
            "$SENTIEON_INSTALL/bin/sentieon util sort -t {threads} -o {output} --sam2bam - &> {log}"


# Sentieon's software uses a pcr indel model to limit false positives
rule indel_modeling:
    input:
        "Assembly/cg50X_SentieonIter0.bam"
    output:
        "Assembly/LinkedReads.pcr.model"
    log:
        "Assembly/LinkedReadsRepeatStat.log"
    threads:
        config['threads']['assembly']
    params:
        seed_lib = config['params']['seed_lib']
    shell:
        "$SENTIEON_INSTALL/bin/sentieon driver -t {threads} "
            "-r {params.seed_lib} -i {input} "
            "--algo RepeatStat {output} &> {log}"


# First assembly step
rule linked_reads_assemble:
    input:
        bam = "Assembly/cg50X_SentieonIter0.bam",
        model = "Assembly/LinkedReads.pcr.model"
    output:
        bam = "Assembly/LinkedReads.hap.bam"
    log:
        "Assembly/LinkedReadsAssemble.log"
    threads:
        config['threads']['assembly']
    params:
        seed_lib = config['params']['seed_lib'],
        read_len = config['params']['read_len']
    shell:
        "$SENTIEON_INSTALL/bin/sentieon driver -t {threads} "
            "-r {params.seed_lib} -i {input.bam} "
            "--algo LinkedReadsAssemble --assemble_read_len {params.read_len} "
            "--pcr_indel_model {input.model} "
            "--coalign_hap 1 --min_map_qual 1 "
            "{output.bam} &> {log}"


# First graph solving step
rule linked_reads_solve:
    input:
        clusters = "Assembly/cg50X_SentieonIter0.bam",
        bam = "Assembly/LinkedReads.hap.bam"
    output:
        graph = "Assembly/LinkedReads.hap.graph",
        reads = "Assembly/LinkedReads.hap.reads",
        filtered = "Assembly/LinkedReads.filtered.reads",
        rescue = "Assembly/LinkedReads.rescue.data"
    log:
        "Assembly/LinkedReadsSolve.log"
    threads:
        config['threads']['assembly']
    params:
        seed_lib = config['params']['seed_lib'],
        range_bc = config['params']['range_barcode']
    shell:
        "$SENTIEON_INSTALL/bin/sentieon driver -t {threads} "
            "-r {params.seed_lib} -i {input.clusters} "
            "--algo LinkedReadsSolve --hap_bam {input.bam} "
            "--hap_graph {output.graph} --reads_in_hap {output.reads} "
            "--reads_filtered {output.filtered} "
            "--rescue_data {output.rescue} "
            "--range_barcode {params.range_bc} &> {log}"


# rescue step for unassembled reads and contigs
rule linked_reads_rescue:
    input:
        bam = "Assembly/cg50X_SentieonIter0.bam",
        rescue = "Assembly/LinkedReads.rescue.data",
        reads = "Assembly/LinkedReads.hap.reads",
        filtered = "Assembly/LinkedReads.filtered.reads"
    output:
        "Assembly/LinkedReads.rescue.bam"
    log:
        "Assembly/LinkedReadsRescue.log"
    threads:
        config['threads']['assembly']
    shell:
        "$SENTIEON_INSTALL/bin/sentieon driver -t {threads} "
            "-i {input.bam} --algo LinkedReadsRescue "
            "--bam_compression 1 --rescue_data {input.rescue} "
            "--reads_in_hap {input.reads} "
            "--reads_filtered {input.filtered} "
            "{output} &> {log}"


# Second assembly step/scaffolding
rule linked_reads_bridge:
    input:
        bam = "Assembly/LinkedReads.rescue.bam",
        data = "Assembly/LinkedReads.rescue.data"
    output:
        "Assembly/LinkedReads.bridge.ret"
    log:
        "Assembly/LinkedReadsBridge.log"
    threads:
        config['threads']['assembly']
    params:
        seed_lib = config['params']['seed_lib']
    shell:
        "$SENTIEON_INSTALL/bin/sentieon driver -t {threads} "
            "-r {params.seed_lib} -i {input.bam} "
            "--algo LinkedReadsBridge --rescue_data {input.data} "
            "{output} &> {log}"


# Second graph solving step
rule linked_reads_solve_2:
    input:
        bam = "Assembly/cg50X_SentieonIter0.bam",
        hap = "Assembly/LinkedReads.hap.bam",
        graph = "Assembly/LinkedReads.hap.graph",
        bridge = "Assembly/LinkedReads.bridge.ret"
    output:
        contigs = expand("Assembly/LinkedReads.contig_{hap}.fasta", hap=[0,1]),
        phase_contigs = expand("Assembly/LinkedReads.contig.phase_{hap}.fasta", hap=[0,1]),
        barcodes = "Assembly/LinkedReads.phase.barcodes"
    log:
        "Assembly/LinkedReadsSolve1.log"
    threads:
        config['threads']['assembly']
    params:
        jaccard = config['params']['min_jaccard'],
        seed_lib = config['params']['seed_lib'],
        range_bc = config['params']['range_barcode']
    shell:
        "$SENTIEON_INSTALL/bin/sentieon driver -t {threads} "
            "-r {params.seed_lib} -i {input.bam} "
            "--algo LinkedReadsSolve --hap_bam {input.hap} "
            "--hap_graph {input.graph} --bridge_ret {input.bridge} "
            "--contig_output Assembly/LinkedReads.contig "
            "--min_scaffold_jaccard {params.jaccard} "
            "--range_barcode {params.range_bc} "
            "--phase_barcodes {output.barcodes} &> {log}"


# run quast on contig assemblies
rule run_quast_contigs:
    input:
        "Assembly/LinkedReads.contig_{hap}.fasta"
    output:
        "Quast_Contigs_{hap}/quast.log"
    threads:
        config['threads']['quast']
    params:
        quast = config['software']['quast'],
        min_contig = config['params']['min_contig']
    shell:
        "{params.quast} {input} --eukaryote "
            "--threads {threads} "
            "--no-snps "
            "--fragmented "
            "--min-contig {params.min_contig} "
            "-o $(dirname {output}) "
            "-s "


# run minimap for each phase hap against truth_haps for binning
rule run_minimap_phased_contigs:
    input:
        fasta = "Assembly/LinkedReads.contig.phase_{hap}.fasta",
        truth_fa = config['samples']['sim_ref_stem'] + "/mergeScaftig_normalized_hap{ref_hap}.fa"
    output:
        "Align_Hap{hap}_Ref{ref_hap}/hap{hap}_ref{ref_hap}.paf.gz"
    params:
        minimap = config['software']['minimap']
    shell:
        "{params.minimap} -x asm5 "
            "--paf-no-hit "
            "--secondary=no "
            "--cs "
            "-y "
            "{input.truth_fa} <(sed 's/=/:i:/g' {input.fasta}) | "
            "gzip -c > {output}"


# get contig N50 using python script
rule contig_phase_n50:
    input:
        "Assembly/LinkedReads.contig.phase_{hap}.fasta"
    output:
        "hap{hap}_phase_n50.txt"
    params:
        phase_script = config['software']['phase_script']
    shell:
        "{params.phase_script} {input} > {output}"


# bin contigs by the phase id (PID), creates temp file of results
rule bin_by_pid:
    input:
        input_fa = expand("Assembly/LinkedReads.contig.phase_{hap}.fasta", hap=[0,1]),
        pafs = expand("Align_Hap{hap}_Ref{ref_hap}/hap{hap}_ref{ref_hap}.paf.gz",
                       hap=[0,1], ref_hap=["A", "B"])
    output:
        "Contig_Binning/pid_binning_results.txt"
    params:
        bin_by_pid = config['software']['bin_by_pid']
    shell:
        "python3 {params.bin_by_pid} "
            "{input.input_fa} "
            "{input.pafs} "
            "{output}"


# maps truth haplotypes to 0 and 1
def get_truth_hap_id(wildcards):
  '''
  Given a truth hap name, get a 0/1 ID
  '''
  ref_hap = wildcards.ref_hap
  if ref_hap == "A":
    return '0'
  elif ref_hap == "B":
    return '1'
  else:
    raise ValueError("Unknown 'truth_hap': " + wildcards.ref_hap)


# returns wildcard for haplotype
def get_contig_hap_id(wildcards):
    '''
    Given a contig hap name, get a 0/1 ID
    '''
    hap = wildcards.hap
    if hap == "0":
        return '0'
    elif hap == "1":
        return '1'
    else:
        raise ValueError("Unknown 'truth_hap': " + wildcards.hap)


# create a file of the appropriate phases for each contig
# creates list of contigs for each original haplotype and truth bin
rule create_subset_files:
    input:
        bin_results = "Contig_Binning/pid_binning_results.txt",
        contig = "Assembly/LinkedReads.contig.phase_{hap}.fasta",
        truth_fa = config['samples']['sim_ref_stem'] + "/mergeScaftig_normalized_hap{ref_hap}.fa"
    output:
        "Contig_Binning/tmp_contigs_hap{hap}_bin{ref_hap}.txt"
    params:
        truth_hap_id = get_truth_hap_id,
        contig_hap_id = get_contig_hap_id
    shell:
        "grep \"^{params.contig_hap_id}[[:space:]]{params.truth_hap_id}\" "
            "{input.bin_results} | "
            "cut -f 3 > {output}"


# subset contigs based on the tmp_contigs_hapN_binZ.txt files
rule subset_bins:
    input:
        bin_results = "Contig_Binning/pid_binning_results.txt",
        zero_contig = "Assembly/LinkedReads.contig.phase_0.fasta",
        zero_subfile = "Contig_Binning/tmp_contigs_hap0_bin{ref_hap}.txt",
        one_contig = "Assembly/LinkedReads.contig.phase_1.fasta",
        one_subfile = "Contig_Binning/tmp_contigs_hap1_bin{ref_hap}.txt"
    output:
        "Contig_Binning/binned_hap{ref_hap}_contigs.fa"
    params:
        seqtk = config['software']['seqtk']
    shell:
        "{params.seqtk} subseq -l 60 "
            "{input.zero_contig} "
            "{input.zero_subfile} >> "
            "{output}; "
        "{params.seqtk} subseq -l 60 "
            "{input.one_contig} "
            "{input.one_subfile} >> "
            "{output}"


# map binned contigs to ref hapA
rule run_minimap_binned_contigs_refA:
    input:
        contigs = "Contig_Binning/binned_hap{ref_hap}_contigs.fa",
        truth_fa = "{}/mergeScaftig_normalized_hapA.fa".format(config['samples']['sim_ref_stem'])
    output:
        paf = "Contig_Binning/Hap{ref_hap}_to_RefA/hap{ref_hap}_refA.paf",
        bam = "Contig_Binning/Hap{ref_hap}_to_RefA/hap{ref_hap}_refA.bam",
        index = "Contig_Binning/Hap{ref_hap}_to_RefA/hap{ref_hap}_refA.bam.bai",
        vcf = "Contig_Binning/Hap{ref_hap}_to_RefA/hap{ref_hap}_refA.vcf"
    params:
        minimap = config['software']['minimap']
    shell:
        "{params.minimap} -x asm5 "
            "--cs "
            "--paf-no-hit "
            "--secondary=no "
            "{input.truth_fa} {input.contigs} > "
            "{output.paf}; "
        "sort -k6,6 -k8,8n {output.paf} | "
            "paftools.js call -l 100 "
            "-L 200 "
            "-f {input.truth_fa} - > {output.vcf}; "
        "{params.minimap} -x asm5 -a "
            "--secondary=no "
            "{input.truth_fa} {input.contigs} | "
            "samtools sort -O bam -o {output.bam}; "
        "samtools index {output.bam}"


# map binned contigs to ref hapB
rule run_minimap_binned_contigs_refB:
    input:
        contigs = "Contig_Binning/binned_hap{ref_hap}_contigs.fa",
        truth_fa = "{}/mergeScaftig_normalized_hapB.fa".format(config['samples']['sim_ref_stem'])
    output:
        paf = "Contig_Binning/Hap{ref_hap}_to_RefB/hap{ref_hap}_refB.paf",
        bam = "Contig_Binning/Hap{ref_hap}_to_RefB/hap{ref_hap}_refB.bam",
        index = "Contig_Binning/Hap{ref_hap}_to_RefB/hap{ref_hap}_refB.bam.bai",
        vcf = "Contig_Binning/Hap{ref_hap}_to_RefB/hap{ref_hap}_refB.vcf"
    params:
        minimap = config['software']['minimap']
    shell:
        "{params.minimap} -x asm5 "
            "--cs "
            "--paf-no-hit "
            "--secondary=no "
            "{input.truth_fa} {input.contigs} > "
            "{output.paf}; "
        "sort -k6,6 -k8,8n {output.paf} | "
            "paftools.js call -l 100 "
            "-L 200 "
            "-f {input.truth_fa} - > {output.vcf}; "
        "{params.minimap} -x asm5 -a "
            "--secondary=no "
            "{input.truth_fa} {input.contigs} | "
            "samtools sort -O bam -o {output.bam}; "
        "samtools index {output.bam}"


# return truth fasta based on wildcards
def get_truth_fa(wildcards):
    return "{}/mergeScaftig_normalized_hap{}.fa".format(config['samples']['sim_ref_stem'], wildcards.ref_hap)


# run quast based on binned contigs
rule run_quast_binned_contigs:
    input:
        contigs = "Contig_Binning/binned_hap{ref_hap}_contigs.fa",
        truth_fa = get_truth_fa
    output:
        "Contig_Binning/Quast_Contigs_Hap{ref_hap}/quast.log"
    threads:
        config['threads']['quast']
    params:
        quast = config['software']['quast'],
        min_contig = config['params']['min_contig']
    shell:
        "{params.quast} {input.contigs} --eukaryote "
            "-R {input.truth_fa} "
            "--threads {threads} "
            "--no-snps "
            "--fragmented "
            "--min-contig {params.min_contig} "
            "-o $(dirname {output}) "
            "-s "


# sort the binned contigs by SID and contig number for scaffolding
# this creates a TSV with a captured contig number and scaffold ID then
# does a numeric sort based on those before transforming the results
# back to fasta format
rule sort_binned_contigs:
    input:
        "Contig_Binning/binned_hap{ref_hap}_contigs.fa"
    output:
        "Contig_Binning/binned_hap{ref_hap}_contigs.sorted.fa"
    params:
        seqkit = config['software']['seqkit'],
        csvtk = config['software']['csvtk'],
        seqtk = config['software']['seqtk']
    shell:
        "cat {input} | "
        "{params.seqkit} fx2tab | "
        "tr \" \" \"\\t\" | "
        "{params.csvtk} -H -t mutate -p \"C([0-9]*)\" | "
        "{params.csvtk} -H -t mutate -p \"SID=([0-9]*)\" | "
        "{params.csvtk} -H -t sort -k 9:n -k 8:n | "
        "awk '{{print \">\"$1\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\n\"$6}}' | "
        "{params.seqtk} seq -l60 - > "
            "{output}"


# scaffold contigs
rule scaffold_contigs:
    input:
        "Contig_Binning/binned_hap{ref_hap}_contigs.sorted.fa"
    output:
        "Scaffolds/scaffolded_hap{ref_hap}_contigs.fa"
    params:
        scaffold = config['software']['scaffold']
    shell:
        "{params.scaffold} {input} > {output}"


# map scaffolded contigs to ref hapA
# creates a VCF, paf and bam file
rule run_minimap_binned_scaffolds_refA:
    input:
        scaffolds = "Scaffolds/scaffolded_hap{ref_hap}_contigs.fa",
        truth_fa = "{}/mergeScaftig_normalized_hapA.fa".format(config['samples']['sim_ref_stem'])
    output:
        paf = "Scaffolds/Hap{ref_hap}_to_RefA/hap{ref_hap}_refA.paf",
        bam = "Scaffolds/Hap{ref_hap}_to_RefA/hap{ref_hap}_refA.bam",
        index = "Scaffolds/Hap{ref_hap}_to_RefA/hap{ref_hap}_refA.bam.bai",
        vcf = "Scaffolds/Hap{ref_hap}_to_RefA/hap{ref_hap}_refA.vcf"
    params:
        minimap = config['software']['minimap']
    shell:
        "{params.minimap} -x asm5 "
            "--cs "
            "--paf-no-hit "
            "--secondary=no "
            "{input.truth_fa} {input.scaffolds} > "
            "{output.paf}; "
        "sort -k6,6 -k8,8n {output.paf} | "
            "paftools.js call -l 100 "
            "-L 200 "
            "-f {input.truth_fa} - > {output.vcf}; "
        "{params.minimap} -x asm5 -a "
            "--secondary=no "
            "{input.truth_fa} {input.scaffolds} | "
            "samtools sort -O bam -o {output.bam}; "
        "samtools index {output.bam}"


# map scaffolded contige to ref hapB
# creates a VCF, paf and bam file
rule run_minimap_binned_scaffolds_refB:
    input:
        scaffolds = "Scaffolds/scaffolded_hap{ref_hap}_contigs.fa",
        truth_fa = "{}/mergeScaftig_normalized_hapB.fa".format(config['samples']['sim_ref_stem'])
    output:
        paf = "Scaffolds/Hap{ref_hap}_to_RefB/hap{ref_hap}_refB.paf",
        bam = "Scaffolds/Hap{ref_hap}_to_RefB/hap{ref_hap}_refB.bam",
        index = "Scaffolds/Hap{ref_hap}_to_RefB/hap{ref_hap}_refB.bam.bai",
        vcf = "Scaffolds/Hap{ref_hap}_to_RefB/hap{ref_hap}_refB.vcf"
    params:
        minimap = config['software']['minimap']
    shell:
        "{params.minimap} -x asm5 "
            "--cs "
            "--paf-no-hit "
            "--secondary=no "
            "{input.truth_fa} {input.scaffolds} > "
            "{output.paf}; "
        "sort -k6,6 -k8,8n {output.paf} | "
            "paftools.js call -l 100 "
            "-L 200 "
            "-f {input.truth_fa} - > {output.vcf}; "
        "{params.minimap} -x asm5 -a "
            "--secondary=no "
            "{input.truth_fa} {input.scaffolds} | "
            "samtools sort -O bam -o {output.bam}; "
        "samtools index {output.bam}"


# run quast with scaffolds against the appropriate truth fasta
rule run_quast_scaffolds:
    input:
        scaffolds = "Scaffolds/scaffolded_hap{ref_hap}_contigs.fa",
        truth_fa = get_truth_fa
    output:
        "Scaffolds/Quast_Scaffolds_Hap{ref_hap}/quast.log"
    threads:
        config['threads']['quast']
    params:
        quast = config['software']['quast'],
        min_contig = config['params']['min_contig']
    shell:
        "{params.quast} {input.scaffolds} --eukaryote "
            "-R {input.truth_fa} "
            "--threads {threads} "
            "--no-snps "
            "--fragmented "
            "--min-contig {params.min_contig} "
            "-o $(dirname {output}) "
            "-s "


# map contigs to a human reference, create a paf file
rule map_contigs_to_hrg:
    input:
        "Contig_Binning/binned_hap{ref_hap}_contigs.fa"
    output:
        "Contig_Binning/SV_Calling/PAF_Alignments/hap{ref_hap}_ref.paf.gz"
    params:
        human_ref = config['params']['human_ref'],
        minimap = config['software']['minimap']
    threads:
        config['threads']['minimap']
    shell:
        "{params.minimap} -cxasm5 "
            "--paf-no-hit "
            "--cs "
            "-r2k "
            "-t {threads} "
            "{params.human_ref} {input} | "
        "gzip -c > {output}"


# call variants in var format and scrape unique alignments
rule get_contig_alignment_bed:
    input:
        "Contig_Binning/SV_Calling/PAF_Alignments/hap{ref_hap}_ref.paf.gz"
    output:
        "Contig_Binning/SV_Calling/PAF_BED_Alignments/hap{ref_hap}_ref.bed"
    shell:
        "zcat {input} | "
        "sort -k6,6 -k8,8n | "
        "paftools.js call - | "
        "grep ^R | "
        "cut -f2- > {output}"


# create diploid coverage bedfile based on contig alignments
rule get_diploid_bed:
    input:
        expand("Contig_Binning/SV_Calling/PAF_BED_Alignments/hap{ref_hap}_ref.bed", ref_hap=['A', 'B'])
    output:
        "Contig_Binning/SV_Calling/PAF_BED_Alignments/diploid_cov.bed"
    params:
        bedtools = config['software']['bedtools']
    shell:
        "{params.bedtools} intersect -a {input[0]} "
            "-b {input[1]} | "
        "grep -v '^chr[0-9]*_' > "
            "{output}"


# map contig to human ref in sam format
rule sam_alignments_to_hrg:
    input:
        "Contig_Binning/binned_hap{ref_hap}_contigs.fa"
    output:
        "Contig_Binning/SV_Calling/SAM_Alignments/hap{ref_hap}_ref.sam.gz"
    params:
        human_ref = config['params']['human_ref'],
        minimap = config['software']['minimap']
    threads:
        config['threads']['minimap']
    shell:
        "{params.minimap} -axasm5 "
            "--paf-no-hit "
            "--cs "
            "-r2k "
            "-t {threads} "
            "{params.human_ref} {input} | "
        "gzip -c > {output}"


# sort and convert the sam to bam
rule sam_to_bam:
    input:
        "Contig_Binning/SV_Calling/SAM_Alignments/hap{ref_hap}_ref.sam.gz"
    output:
        "Contig_Binning/SV_Calling/BAM_Alignments/hap{ref_hap}_ref.bam"
    params:
        sam_flt = config['software']['sam_flt']
    shell:
        "{params.sam_flt} {input} | "
        "samtools sort -m4G "
            "-@4 "
            "-o {output}"


# create a paired VCF file
rule pair_vcfs:
    input:
        expand("Contig_Binning/SV_Calling/BAM_Alignments/hap{ref_hap}_ref.bam", ref_hap = ['A', 'B'])
    output:
        "Contig_Binning/SV_Calling/Paired_VCF/query_contigs_pair.vcf.gz"
    params:
        htsbox = config['software']['htsbox'],
        human_ref = config['params']['human_ref']
    shell:
        "{params.htsbox} pileup "
            "-q5 "
            "-evcf "
            "{params.human_ref} "
            "{input} | "
        "{params.htsbox} bgzip > {output}"


# phase vcfs with vcf_pair
rule phase_vcfs:
    input:
        "Contig_Binning/SV_Calling/Paired_VCF/query_contigs_pair.vcf.gz"
    output:
        "Contig_Binning/SV_Calling/Phased_VCFs/query_contigs_phased.vcf.gz"
    params:
        vcf_pair = config['software']['vcf_pair'],
        htsbox = config['software']['htsbox']
    shell:
        "{params.vcf_pair} {input} | "
        "awk '{{ if ($0 ~ /^#/ || ($4 ~ /^[ACGT]/ && $5 ~ /^[ACGT]/ )) print }}' | "
        "{params.htsbox} bgzip > {output}"


# index the phased VCF
rule index_phased_vcf:
    input:
        "Contig_Binning/SV_Calling/Phased_VCFs/query_contigs_phased.vcf.gz"
    output:
        "Contig_Binning/SV_Calling/Phased_VCFs/query_contigs_phased.vcf.gz.tbi"
    shell:
        "tabix -p vcf {input}"


# compare the VCF results to benchmarked SVs made from the truth haplotypes
# uses truvari for the comparison
rule benchmark_svs:
    input:
        vcf = "Contig_Binning/SV_Calling/Phased_VCFs/query_contigs_phased.vcf.gz",
        index = "Contig_Binning/SV_Calling/Phased_VCFs/query_contigs_phased.vcf.gz.tbi"
    output:
        "Contig_Binning/SV_Calling/SV_Comparison/summary.txt"
    params:
        truth_vcf = config['params']['sv_truth'],
        ref = config['params']['human_ref']
    shell:
        "source ~/Virtual_Envs/Truvari/bin/activate; "
        "rm -r $(dirname {output}); "
        "truvari -b {params.truth_vcf} "
            "-c {input.vcf} "
            "-o $( dirname {output} ) "
            "-f {params.ref}"


# benchmark again using only diploid coverage sites
rule benchmark_svs_diploid_cov:
    input:
        vcf = "Contig_Binning/SV_Calling/Phased_VCFs/query_contigs_phased.vcf.gz",
        index = "Contig_Binning/SV_Calling/Phased_VCFs/query_contigs_phased.vcf.gz.tbi",
        bed = "Contig_Binning/SV_Calling/PAF_BED_Alignments/diploid_cov.bed"
    output:
        "Contig_Binning/SV_Calling/SV_Comparison_Diploid/summary.txt"
    params:
        truth_vcf = config['params']['sv_truth'],
        ref = config['params']['human_ref']
    shell:
        "source ~/Virtual_Envs/Truvari/bin/activate; "
        "rm -r $(dirname {output}); "
        "truvari -b {params.truth_vcf} "
            "-c {input.vcf} "
            "-o $( dirname {output} ) "
            "--includebed {input.bed} "
            "-f {params.ref}"


# stratify the SVs using stratification_test.py
rule stratify_svs:
    input:
        "Contig_Binning/SV_Calling/SV_Comparison/summary.txt"
    output:
        "Contig_Binning/SV_Calling/Stratification_Results/results.png"
    params:
        strat_tsv = config['params']['strat_tsv'],
        strat_bed_dir = config['params']['strat_bed_dir'],
        strat_svs = config['software']['strat_svs_script']
    shell:
        "{params.strat_svs} -t {params.strat_tsv} "
            "-b {params.strat_bed_dir} "
            "-o $( dirname {output} ) "
            "$(dirname {input} )"
