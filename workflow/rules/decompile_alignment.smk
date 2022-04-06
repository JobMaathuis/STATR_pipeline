rule convert_sam_to_bam:
    message: 'converting .sam to .bam of {wildcards.sample}'
    input:
        config['wdir'] + config['result'] + '2.aligned/{sample}.sam'
    benchmark:
        'benchmarks/convert_sam_to_bam_{sample}.txt'
    output:
        config['wdir'] + config['result'] + '3.decompiled/{sample}.bam'
    shell:
        'mkdir -p ' + config['wdir'] + config['result'] + '3.decompiled;'
        'samtools view -bS {input} > {output}'


rule sort_bam:
    message: 'sorting .bam file of {wildcards.sample}'
    input:
        config['wdir'] + config['result'] + '3.decompiled/{sample}.bam'
    benchmark:
        'benchmarks/sort_bam_{sample}.txt'
    output:
        config['wdir'] + config['result'] + '3.decompiled/{sample}_sorted.bam'
    shell:
        'samtools sort {input} -o {output}'


rule convert_bam_to_bed:
    message: 'converting .bam to .bed of {wildcards.sample}'
    input:
        config['wdir'] + config['result'] + '3.decompiled/{sample}_sorted.bam'
    benchmark:
        'benchmarks/convert_bam_to_bed_{sample}.txt'
    output:
        config['wdir'] + config['result'] + '3.decompiled/{sample}.bed'
    shell:
        'bedtools bamtobed -i {input} > {output}'