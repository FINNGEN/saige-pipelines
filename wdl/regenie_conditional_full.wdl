version 1.0

workflow conditional_analysis {

  input {
    String docker
    File phenos_to_cond
   
    String chr_col
    String pos_col
    String ref_col
    String alt_col
    String pval_col
    Float conditioning_pval_threshold
    Float locus_pval_threshold  
    
    Boolean add_chr
  }
  
  Array[Array[String]] pheno_data = read_tsv(phenos_to_cond)
  # go through phenos
  scatter (p in pheno_data) {
    #get hits under pval threshold
    call extract_cond_regions {
      input: pheno=p[0], summaryfile=p[1], pval_col=pval_col, pval_threshold=locus_pval_threshold, docker=docker
    }
  }
}

task extract_cond_regions {

  input {
    File summaryfile
    String pheno
    String pval_col
    String docker
    Float pval_threshold
  }
  
  String outfile=basename(summaryfile,".gz" )+ "-gw_sig.tsv"
  
  command <<<
    zcat ~{summaryfile} | awk -v pcol=~{pval_col} -v pheno=~{pheno} -v pth=~{pval_threshold}  \
    'BEGIN{ FS="\t"; OFS="\t"; p_idx=0} NR==1{ for(i=1;i<=NF; i++) { if($i==pcol) p_idx=i  }; if(p_idx==0) {print "Given p-val column not found" > "/dev/stderr"; exit 1} print "PHENO",$0     } \
    NR>1{ if( $(p_idx)<pth && $(p_idx)>0) print pheno,$0} ' > ~{outfile}
    >>>
    
    output {
      File gw_sig_res = outfile
    }
    
    runtime {
      cpu: "1"
      docker: "${docker}"
      memory: "2 GB"
      disks: "local-disk 10 HDD"
      zones: "europe-west1-b europe-west1-c europe-west1-d"
      preemptible: "1"
    }
  }
  




  
