rule index_genome:
    message: 'indexing ' + config['genome']
    input:
        genome = config['wdir'] + config['files'] + config['genome'] + '.fa'
    output:
        temp(touch('index.done'))
    benchmark:
        'benchmarks/index_genome.txt'
    params:
        index = config['wdir'] + config['result'] + '2.aligned/' + config['genome']
    shell:
        'mkdir -p ' + config['wdir'] + config['result'] + '2.aligned;'
        'bowtie2-build -f {input.genome} {params.index}'


rule align_reads:
    message: 'aliging reads of {wildcards.sample}'
    input: 
        check = 'index.done',
        trimmed = config['wdir'] + config['result'] + '1.trimmed/{sample}_trimmed' + config['fastq']
    output:
        config['wdir'] + config['result'] + '2.aligned/{sample}.sam'
    benchmark:
        'benchmarks/align_reads_{sample}.txt'
    threads: 4
    params:
        index = config['wdir'] + config['result'] + '2.aligned/' + config['genome']
    shell:
        'bowtie2 -q -U {input.trimmed} -x {params.index} -S {output} --local --threads {threads}'