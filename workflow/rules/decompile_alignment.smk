rule convert_sam_to_bam:
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
    message: 'converting .bam to .bed of {wildcards.sample}'
    input:
        config['wdir'] + config['result'] + '3.decompiled/{sample}_sorted.bam'
    log:
        config['wdir'] + 'logs/3.decompiled/{sample}_bam_to_bed.log'
    output:
        config['wdir'] + config['result'] + '3.decompiled/{sample}.bed'
    shell:
        '(bedtools bamtobed -i {input} > {output}) 2> {log}'