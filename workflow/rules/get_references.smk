rule get_genome:
    output:
        config["genome"]["outdir"]+"/"+config["genome"]["species"]+"/genome.fasta"
    params:
        species=config["genome"]["species"],
        build=config["genome"]["build"],
        release=config["genome"]["release"],
        datatype=config["genome"]["datatype"]
    log:
        "logs/get_genome.log"
    wrapper:
        "v1.31.1/bio/reference/ensembl-sequence"

rule get_gtf:
    output:
        config["genome"]["outdir"]+"/"+config["genome"]["species"]+"/genome.gtf"
    params:
        species=config["genome"]["species"],
        build=config["genome"]["build"],
        release=config["genome"]["release"],
        flavor=""
    log:
        "logs/get_gtf.log"
    wrapper:
        "v1.31.1/bio/reference/ensembl-annotation"