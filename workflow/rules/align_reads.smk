rule index_genome:
    """ Indexes the genome using Bowtie2 (which uses FM index) """
    message: 'indexing ' + config['genome']
    input:
        genome = config['wdir'] + config['files'] + config['genome'] + '.fa'
    output:
        touch('index.done')
    log:
        config['wdir'] + 'logs/2.aligned/index_genome.log'
    params:
        index = config['wdir'] + config['result'] + '2.aligned/' + config['genome']
    shell:
        'mkdir -p ' + config['wdir'] + config['result'] + '2.aligned;'
        '(bowtie2-build -f {input.genome} {params.index}) 2> {log}'


rule align_reads:
    """ Aligns the trimmed reads to the indexed genome, which creates a .sam file for every sample """
    message: 'aliging reads of {wildcards.sample}'
    input: 
        check = 'index.done',
        trimmed = config['wdir'] + config['result'] + '1.trimmed/{sample}_trimmed' + config['fastq']
    output:
        config['wdir'] + config['result'] + '2.aligned/{sample}.sam'
    log:
        config['wdir'] + 'logs/2.aligned/{sample}_aligning.log'
    threads: 4
    params:
        index = config['wdir'] + config['result'] + '2.aligned/' + config['genome']
    shell:
        '(bowtie2 -q -U {input.trimmed} -x {params.index} -S {output} --local --threads {threads}) 2> {log}'