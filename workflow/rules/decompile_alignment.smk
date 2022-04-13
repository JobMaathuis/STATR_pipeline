rule convert_sam_to_bam:
    """ Converts the .sam files to .bam files, which is a binary format and advantageous for computer programs """
    message: 'converting .sam to .bam of {wildcards.sample}'
    input:
        config['wdir'] + config['result'] + '2.aligned/{sample}.sam'
    log:
        config['wdir'] + 'logs/3.decompiled/{sample}_sam_to_bam.log'
    output:
        config['wdir'] + config['result'] + '3.decompiled/{sample}.bam'
    shell:
        'mkdir -p ' + config['wdir'] + config['result'] + '3.decompiled;'
        '(samtools view -bS {input} > {output}) 2> {log}'


rule sort_bam:
    """ Sorts the .bam files by leftmost coordinates using samtools """
    message: 'sorting .bam file of {wildcards.sample}'
    input:
        config['wdir'] + config['result'] + '3.decompiled/{sample}.bam'
    log:
        config['wdir'] + 'logs/3.decompiled/{sample}_sort_bam.log'
    output:
        config['wdir'] + config['result'] + '3.decompiled/{sample}_sorted.bam'
    shell:
        '(samtools sort {input} -o {output}) 2> {log}'


rule convert_bam_to_bed:
    """ In order to store the genomic regions by coordinates and annotation the .bam files are converted to .bed files using bedtools """
    message: 'converting .bam to .bed of {wildcards.sample}'
    input:
        config['wdir'] + config['result'] + '3.decompiled/{sample}_sorted.bam'
    log:
        config['wdir'] + 'logs/3.decompiled/{sample}_bam_to_bed.log'
    output:
        config['wdir'] + config['result'] + '3.decompiled/{sample}.bed'
    shell:
        '(bedtools bamtobed -i {input} > {output}) 2> {log}'