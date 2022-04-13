rule trim_reads:
    """ Trims the input files where it removes the adapter and other illumina-specific sequences """
    message: 'trimming reads of {wildcards.sample}'
    input:
        config['wdir'] + config['files'] + '{sample}' + config['fastq']
    log:
        config['wdir'] + 'logs/1.trimmed/{sample}_trimming.log'
    output:
        config['wdir'] + config['result'] + '1.trimmed/{sample}_trimmed' + config['fastq']
    threads: 2
    shell:
        'mkdir -p ' + config['wdir'] + config['result'] + '1.trimmed;'
        '(java -jar ' + config['wdir'] + config['resources'] + config['trimmomatic'] + ' SE -threads {threads} -phred33 {input} {output} '
        'ILLUMINACLIP:' + config['wdir'] + config['files'] + config['adapter'] + ':2:30:6 SLIDINGWINDOW:4:15 MINLEN:12) 2> {log}'