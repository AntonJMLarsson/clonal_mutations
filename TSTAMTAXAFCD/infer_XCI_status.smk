configfile: "config.yaml"
RUNS, = glob_wildcards("cellranger_files/{run}/sample_alignments.bam")

rule all: 
    input: "vireo_out/GT_donors.vireo.vcf.gz"

rule make_barcode_list:
    input: "cellranger_files/{run}/sample_filtered_barcodes.csv"
    output: "cellranger_files/{run}/barcode_list.txt"
    shell: """awk -F ',' '{{print $2}}' {input} > {output}"""

rule modify_barcode_list:
    input: "cellranger_files/{run}/barcode_list.txt"
    output: "cellranger_files/{run}/barcode_list.modified.txt"
    shell: """awk '{{print $0"_{wildcards.run}"}}' {input} > {output}"""

rule cat_barcode_lists:
    input: expand("cellranger_files/{run}/barcode_list.modified.txt", run=RUNS)
    output: "{}_full_barcode_list.txt".format(config["sample"])
    shell: "cat {input} > {output}"
rule extract_reads_X:
    input: bam = "cellranger_files/{run}/sample_alignments.bam", barcodes = "cellranger_files/{run}/barcode_list.txt"
    output: temp("cellranger_files/{run}/sample_alignments.X.bam")
    shell: "python3 {config[src_dir]}/extract_reads_X.py -i {input.bam} -o {output} --contig {config[contig_X]} --cells {input.barcodes} --exon-tag {config[exon_tag]} --intron-tag {config[intron_tag]} --barcode-tag {config[barcode_tag]} --suffix {wildcards.run}"

rule concatenate_bams:
    input: expand("cellranger_files/{run}/sample_alignments.X.bam", run = RUNS)
    output: "{}.{}.bam".format(config["sample"], config["contig_X"])
    shell: "samtools cat -o {output} {input}"

rule sort_bam_X:
    input:
        "{}.{}.bam".format(config["sample"], config["contig_X"])
    output:
        "{}.{}.sorted.bam".format(config["sample"], config["contig_X"])
    threads: config["threads"]
    shell: "samtools sort -@ {threads} -m 2G -o {output} {input}"

rule index_bam_X:
    input: "{}.{}.sorted.bam".format(config["sample"], config["contig_X"])
    output: "{}.{}.sorted.bam.bai".format(config["sample"], config["contig_X"])
    shell: "samtools index {input}"

rule cellsnp:
    input: bam = "{}.{}.sorted.bam".format(config["sample"], config["contig_X"]), bai = "{}.{}.sorted.bam.bai".format(config["sample"], config["contig_X"]), barcodes = "{}_full_barcode_list.txt".format(config["sample"])
    output: "cellsnp/cellSNP.base.vcf.gz"
    threads: config["threads"]
    shell: "cellsnp-lite -s {input.bam} -b {input.barcodes} -O cellsnp -R {config[SNP_file]} -p {threads} --minCOUNT 5 --minMAF {config[MAF]} --cellTAG {config[barcode_tag]} --gzip --UMItag None"

rule vireo:
    input: "cellsnp/cellSNP.base.vcf.gz"
    output: "vireo_out/GT_donors.vireo.vcf.gz"
    threads: config["threads"]
    shell: "vireo -c cellsnp/ -N 2 -o vireo_out --noDoublet -p {threads}"
