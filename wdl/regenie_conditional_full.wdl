version 1.0

workflow conditional_analysis {

  input {
    String docker
    File phenos_to_cond

    String mlogp_col
    Float conditioning_mlogp_threshold
    Float locus_mlogp_threshold  
    Array[String] chroms 
  }
  
  Array[String] pheno_data = read_lines(phenos_to_cond)
  # go through phenos
  scatter (p in pheno_data) {
    #get hits under pval threshold
    call extract_cond_regions {
      input: pheno=p, mlogp_threshold = locus_mlogp_threshold, docker=docker,mlogp_col = mlogp_col,chroms=chroms
    }
  }
  call merge_regions {input: docker =docker,hits=extract_cond_regions.gw_sig_res}
}

  
task extract_cond_regions {
  
  input {
    
    String pheno
    String summary_root
    String region_root
    String chr_col
    String pos_col
    String ref_col
    String alt_col
    String mlogp_col
    Array[String] chroms 

    String docker
    Float mlogp_threshold

  }

  File region = sub(region_root,"PHENO",pheno)
  File sumstats = sub(summary_root,"PHENO",pheno)
  File index = sumstats + ".tbi"
  Int disk_size = ceil(size(sumstats,'GB')) + ceil(size(region,'GB')) + 1
  
  String outfile= pheno + "_sig_hits.txt" 
  
  command <<<

    python3 /scripts/filter_hits_regions.py --sumstats ~{sumstats} --regions ~{region} \
    --pheno ~{pheno} --pval_threshold ~{mlogp_threshold} \
    --pos_col ~{pos_col} --chr_col ~{chr_col} --ref_col ~{ref_col} --alt_col ~{alt_col} --mlogp_col ~{mlogp_col}  --chroms ~{sep=" " chroms} --out ./ --log info
  >>>
  
  output {
    File gw_sig_res = outfile
  }
  
  runtime {
    cpu: "1"
    docker: "${docker}"
    memory: "2 GB"
    disks: "local-disk ${disk_size} HDD"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: "1"
  }
}


task merge_regions {

  input {
    Array[File] hits
    String docker
  }

  String outfile = "regions.txt"
  command <<<
    while read f; do cat $f >> ~{outfile}; done <~{write_lines(hits)}
  >>>

  output {
    File regions = outfile
  }
  
  runtime {
    cpu: "1"
    docker: "${docker}"
    memory: "2 GB"
    disks: "local-disk 2 HDD"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: "1"
  }
}
