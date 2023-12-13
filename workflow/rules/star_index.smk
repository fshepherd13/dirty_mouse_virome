rule star_index:
    input:
        fasta = rules.get_genome.output,
        gtf = rules.get_gtf.output
    output:
        directory(config["genome"]["outdir"]+"/"+config["genome"]["species"]+"/index")
    params:
        mem=40000000000,
        species=config["genome"]["species"],
        build=config["genome"]["build"],
        release=config["genome"]["release"],
        datatype=config["genome"]["datatype"]
    message:
        "Testing STAR index"
    threads: 16
    conda:
        "../envs/star.yml"
    log:
        "logs/star_index.log"
    shell:
        """
        mkdir -p {output}
        STAR \
            --runMode genomeGenerate \
            --runThreadN {threads} \
            --limitGenomeGenerateRAM {params.mem} \
            --genomeDir {output} \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang 150
        
        echo 'Creating genome index for species {params.species}, build {params.build}, release {params.release}, datatype {params.datatype}' > {output}/readme.log
        """