rule filter_contigs:
    input:
        expand("../results/trinity/{exp}_unmapped_denovo_contigs.fasta", exp=EXP)
    output:
        expand("../results/trinity/{exp}_filtered_contigs.fasta", exp=EXP)
    conda:
        "../envs/bbmap.yml"
    shell:
        '''
        reformat.sh in={input} out={output} minlength=500
        '''