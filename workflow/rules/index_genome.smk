rule download_genome:
    message: 'downloading genome'
    input:
        coni

rule index_genome:
    message: 'indexing genome´
    input:
        genome = config['wdir'] + config['files'] + 'genome' + '.fa'
    output:
        temp(touch('index.done'))
    benchmark:
        'benchmarks/index_genome.txt'
    params:
        index = config['wdir'] + config['result'] + '2.aligned/' + 'genome'
    shell:
        'mkdir -p ' + config['wdir'] + config['result'] + '2.aligned;'
        'bowtie2-build -f {input.genome} {params.index}'