rule summarize_transcripts:
    input:
        lineages = rules.assign_lineages.output,
        salmon_files = expand("../results/salmon_quant/{sample}/quant.sf", sample=SAMPLES)
    output:
        species_counts = expand("../results/final/{exp}_sp_abundance.csv", exp=EXP),
        family_counts = expand("../results/final/{exp}_family_abundance.csv", exp=EXP)
    log:
        "logs/summarise_transcripts.log"
    conda:
        "../envs/summary.yml"
    script:
        "../scripts/summarise_transcripts.R"