rule compute_coverage:
    """ Computes the coverage of each of the given genome features """ 
    message: 'computing coverage of {wildcards.sample}'
    input:
        bed = config['wdir'] + config['result'] + '3.decompiled/{sample}.bed',
        gff = config['wdir'] + config['result'] + '4.meta_analysis/' + config['genome'] + '.gff'
    log:
        config['wdir'] + 'logs/5.differential_expression/{sample}_coverage.log'
    output:
        config['wdir'] + config['result'] + '5.differential_expression/{sample}.cov'
    shell:
        'mkdir -p ' + config['wdir'] + config['result'] + '5.differential_expression;'
        '(bedtools coverage -a {input.gff} -b {input.bed} -s > {output}) 2> {log}'


rule format_deseq_input:
    """ Combines all of the samples in one file in a format that can be used by DESeq2 """
    message: 'formating deseq inputs'
    input:
        expand(config['wdir'] + config['result'] + '5.differential_expression/{sample}.cov', sample=config['samples'])
    log:
        config['wdir'] + 'logs/5.differential_expression/format_deseq_input.log'
    output:
        config['wdir'] + config['result'] + '5.differential_expression/DESeq_input.txt'
    shell:
        '(python3 ' + config['wdir'] + config['resources'] + config['py_scripts'] + 'FormatDESeqInput.py -i {input} -o {output}) 2> {log}'


rule run_deseq:
    """ Runs a DESeq2 analysis on the samples, which returns a PCA plot and a heatmap as final output """
    message: "running DESeq"
    input:
        deseq = config['wdir'] + config['result'] + '5.differential_expression/DESeq_input.txt',
    log:
        config['wdir'] + 'logs/5.differential_expression/run_deseq.log'
    output:
        touch('deseq.done')
    params:
        design = config['wdir'] + config['files'] + 'Design_sheet.txt'
    shell:
        '(Rscript ' + config['wdir'] + config['resources'] + config['r_scripts'] + 'RunDESeq.R -i {input.deseq} -d {params.design} -o ' + config['wdir'] + config['result'] + ') 2> {log}'


rule move_txt_files:
    """ Moves the less interesting .txt files to its corresponding directory """
    message: 'moving .txt files'
    input: 
        'deseq.done'
    log:
        config['wdir'] + 'logs/5.differential_expression/move_txt_files.log'
    output:
        touch('all.done')
    shell:
        '(mv ' + config['wdir'] + config['result'] + '*.txt ' + config['wdir'] + config['result'] + '5.differential_expression/) 2> {log}'
