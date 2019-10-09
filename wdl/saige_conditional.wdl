
task condition {
    File null_file
    File variance_ratio_file
    File bgen_file
    File bgen_index
    File sample_file
    String chr
    Int start_pos
    Int end_pos
    String condition_on
    Float p_threshold

    String docker
    Int cpu
    Int mem
    Int disk

    command <<<
        do_condition_locus.py ${null_file} ${variance_ratio_file} \
            ${bgen_file} ${bgen_index} ${sample_file} ${chr} ${start_pos} \
            ${end_pos} ${condition_on} --p_threshold ${p_threshold} --saige_script "Rscript /usr/local/bin/step2_SPAtests.R"
    >>>

    output {
        Array[File] cond_runs = glob("*.conditional")
        Array[File] summary = glob("*independent.snps")
        Array[File] log = glob("*.conditional.log")
    }

    runtime {
        docker: "${docker}"
        cpu: "${cpu}"
        memory: "16 GB"
        disks: "local-disk ${disk} HDD"
        zones: "europe-west1-b"
        preemptible: 2
        noAddress: true
        maxRetries: 0
    }
}

workflow conditional_analysis {
    File loci_to_cond_conf

    String docker
    Int cpu=1
    Int mem=128
    Int disk=20

    Array[Array[String]] loci_to_cond = read_tsv(loci_to_cond_conf)
    scatter (locus in loci_to_cond) {
        call condition {
            input: null_file=locus[0], variance_ratio_file=locus[1],
                bgen_file=locus[2], bgen_index=locus[3], sample_file=locus[4],
                chr=locus[5], start_pos=locus[6], end_pos=locus[7],
                condition_on=locus[8], p_threshold=locus[9], docker=docker,
                cpu=cpu, mem=mem, disk=disk
        }
    }

}
