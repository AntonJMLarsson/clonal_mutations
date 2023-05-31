configfile: "config.yaml"

rule all: 
    input: "vireo_out/GT_donors.vireo.vcf.gz"


rule extract_reads_X:
    input: 
        config["bamfile"]
    output:
        temp("{}.{}.bam".format(config["sample"], config["contig_X"]))
    shell: "python3 {config[src_dir]}/extract_reads_X.py -i {input} -o {output} --contig {config[contig_X]} --cells {config[barcodes]}"

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
    input: bam = "{}.{}.sorted.bam".format(config["sample"], config["contig_X"]), bai = "{}.{}.sorted.bam.bai".format(config["sample"], config["contig_X"])
    output: "cellsnp/cellSNP.base.vcf.gz"
    threads: config["threads"]
    shell: "cellsnp-lite -s {input.bam} -b {config[barcodes]} -O cellsnp -R {config[SNP_file]} -p {threads} --minCOUNT 5 --minMAF {config[MAF]} --cellTAG BC --gzip --UMItag None"

rule vireo:
    input: "cellsnp/cellSNP.base.vcf.gz"
    output: "vireo_out/GT_donors.vireo.vcf.gz"
    threads: config["threads"]
    shell: "vireo -c cellsnp/ -N 2 -o vireo_out --noDoublet -p {threads}"