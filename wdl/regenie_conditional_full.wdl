version 1.0

workflow conditional_analysis {

  input {
    String docker
    File phenos_to_cond

    String pval_col
    Float conditioning_pval_threshold
    Float locus_pval_threshold  
    
  }
  
  Array[String] pheno_data = read_lines(phenos_to_cond)
  # go through phenos
  scatter (p in pheno_data) {
    #get hits under pval threshold
    call extract_cond_regions {
      input: pheno=p, pval_threshold=locus_pval_threshold, docker=docker,pval_col=pval_col
    }
  }
  call gather_hits {
    input: docker = docker,sig_vars = extract_cond_regions.gw_sig_res,p_val_threshold = locus_pval_threshold,conditioning_pval_threshold = conditioning_pval_threshold,pval_col = pval_col}
}

task gather_hits {

  input {
    Array[File] sig_vars

    String pval_col
    String pos_col
    String ref_col
    String alt_col
    String chr_col
    String docker
    Float p_val_threshold
    Float conditioning_pval_threshold
    Boolean add_chr
    String excludes=""

  }

  Int disk_size = ceil(size(sig_vars,'GB'))*2 + 2
  String excl_switch =  if excludes!="" then "--exclude_regions \"" + excludes + "\"" else ""
  String out_file = "regenie_hits.conf"

  command <<<
    cat <( head -n 1 ~{sig_vars[0]}) <( tail -q -n+2 ~{sep=' ' sig_vars} ) > all_vars
    python3 /scripts/generate_conditional_analysis_config.py --p_threshold ~{p_val_threshold} --p_condition_threshold ~{conditioning_pval_threshold}  ~{out_file} vars all_vars \
    --pos_col ~{pos_col} --chr_col "~{chr_col}" --ref_col ~{ref_col} --alt_col ~{alt_col} --p_col ~{pval_col} ~{excl_switch} \
    ~{true='--add_chr' false=' ' add_chr}    
    >>>

  output {
    File all_hits = out_file + ".merged"
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
  
task extract_cond_regions {
  
  input {
    
    String pheno
    String summary_root
    String pval_col

    String docker
    Float pval_threshold
    File? rsid_filter  

  }

  File summaryfile = sub(summary_root,"PHENO",pheno)
  String outfile=basename(summaryfile,".gz" )+ "-gw_sig.tsv"
  
  command <<<
    echo ~{summary_root} ~{pheno}
    zcat ~{summaryfile} | awk -v pcol=~{pval_col} -v pheno=~{pheno} -v pth=~{pval_threshold}  \
    'BEGIN{ FS="\t"; OFS="\t"; p_idx=0} NR==1{ for(i=1;i<=NF; i++) { if($i==pcol) p_idx=i  }; if(p_idx==0) {print "Given p-val column not found" > "/dev/stderr"; exit 1} print "PHENO",$0     } \
    NR>1{ if( $(p_idx)<pth && $(p_idx)>0) print pheno,$0} ' ~{if defined(rsid_filter) then " | { head -n1 ; grep -wf ~{rsid_filter}||true;}" else " "} > ~{outfile}
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





  
