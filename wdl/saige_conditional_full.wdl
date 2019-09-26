import "saige_conditional.wdl" as cond

task extract_cond_regions {
    File summaryfile
    String pheno
    String pval_col
    String docker
    Float pval_threshold
    String outfile=basename(summaryfile,".gz" )+ "-gw_sig.tsv"

    command <<<
        zcat ${summaryfile} | awk -v pcol=${pval_col} -v pheno=${pheno} -v pth=${pval_threshold}  \
        'BEGIN{ FS="\t"; OFS="\t"; p_idx=0} NR==1{ for(i=1;i<=NF; i++) { if($i==pcol) p_idx=i  }; if(p_idx==0) {print "Given p-val column not found" > "/dev/stderr"; exit 1} print "PHENO",$0     } \
            NR>1{ if( $(p_idx)<pth && $(p_idx)>0) print pheno,$0} ' > ${outfile}
    >>>

    output {
        File gw_sig_res = outfile
    }

    runtime {
        cpu: "1"
        docker: "${docker}"
        memory: "2 GB"
        disks: "local-disk 10 HDD"
        zones: "europe-west1-b"
        preemptible: "1"
    }

}


task generate_conditional_analysis_config {
    Array[File] sig_vars
    String null_file_mask
    String var_ratio_mask
    String bgen_file_mask
    String sample_file
    String output_conf
    String docker
    Int padding
    String pval_col
    String chr_col
    String pos_col
    String ref_col
    String alt_col
    Float p_val_threshold
    Float conditioning_pval_threshold

    String excludes=""

    String excl_switch =  if excludes!="" then "--exclude_regions " + excludes else ""

    command <<<
        cat <( head -n 1 ${sig_vars[0]}) <( tail -q -n+2 ${sep=' ' sig_vars} ) > all_vars
        generate_conditional_analysis_config.py --p_threshold ${p_val_threshold} --p_condition_threshold ${conditioning_pval_threshold} \
         ${null_file_mask} ${var_ratio_mask} ${bgen_file_mask} ${sample_file} ${output_conf} vars all_vars --locus_padding ${padding} \
         --pos_col ${pos_col} --chr_col ${chr_col} --ref_col ${ref_col} --alt_col ${alt_col} --p_col ${pval_col} ${excl_switch}
    >>>

    output {
        File conf="${output_conf}.merged"
        File variants=output_conf
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: "2 GB"
        disks: "local-disk 10 HDD"
        zones: "europe-west1-b"
        preemptible: "1"
    }
}

workflow conditional_analysis {
    File phenos_to_cond
    Float locus_pval_threshold
    Float conditioning_pval_threshold
    String chr_col
    String pos_col
    String ref_col
    String alt_col
    String pval_col
    String docker

    Array[Array[String]] pheno_to_cond = read_tsv(phenos_to_cond)
    scatter (p in pheno_to_cond) {
        call extract_cond_regions {
            input: pheno=p[0], summaryfile=p[1], pval_col=pval_col, pval_threshold=locus_pval_threshold, docker=docker
        }
    }

    call generate_conditional_analysis_config {
        input: sig_vars=extract_cond_regions.gw_sig_res, docker=docker,pval_col=pval_col,
            p_val_threshold=locus_pval_threshold, conditioning_pval_threshold=conditioning_pval_threshold,
            chr_col=chr_col,pos_col=pos_col, ref_col=ref_col, alt_col=alt_col
    }

    call cond.conditional_analysis {
        input: loci_to_cond_conf=generate_conditional_analysis_config.conf, docker=docker,
    }

}
