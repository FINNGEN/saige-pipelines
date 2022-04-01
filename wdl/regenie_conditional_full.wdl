version 1.0

workflow conditional_analysis {

  input {
    String docker
    File phenos_to_cond

    String mlogp_col
    Float conditioning_mlogp_threshold
    Float locus_mlogp_threshold  
    Array[String] chroms
    Array[String] covariates
    File pheno_file

  }

  # returns covariate string for each pheno
  call filter_covariates {input: docker=docker,pheno_file=pheno_file,pheno_list = phenos_to_cond,covariates=covariates}  
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


task filter_covariates {

  input {
    File pheno_file
    Array[String] covariates
    File pheno_list
    String docker
    Int threshold_cov_count
    }

    String outfile = "./pheno_cov_map_" + threshold_cov_count + ".txt"
    Int disk_size = ceil(size(pheno_file,'GB')) + 2
    
    command <<<

      python3 <<CODE

      import pandas as pd
      import numpy as np

      #read in phenos as list of phenos regardless
      tot_phenos = []
      phenos_groups = []
      with open('~{pheno_list}') as i:
          for line in i:
              phenos = line.strip().split()
              phenos_groups.append(phenos)
              tot_phenos += phenos    

      #read in phenos mapping all valid entries to 1 and NAs to 0
      pheno_df= pd.read_csv('~{pheno_file}',sep='\t',usecols=tot_phenos).notna().astype(int)
      print(pheno_df)
      # read in covariates getting absolute values (handles PCs)
      covariates= '~{sep="," covariates}'.split(',')
      cov_df= pd.read_csv('~{pheno_file}',sep='\t',usecols=covariates).abs()
      print(cov_df)

      # now for each pheno calculate product of each covariate with itself
      
      with open('~{outfile}','wt') as o,open('~{outfile}'.replace('.txt','.err.txt'),'wt') as tmp_err:
          for i,pheno_list in enumerate(phenos_groups):
              pheno_name = ','.join(pheno_list)
              # for each group of phenos (possibly a single one) multiply all covs and pheno column and count how many non 0 entries are there: this means that the entry has a valid pheno and a non null covariates.
              df = pd.DataFrame()
              for pheno in pheno_list:
                  m = (cov_df.mul(pheno_df[pheno],0)>0).sum().to_frame(pheno)
                  df = pd.concat([df,m],axis =1)

              print(f"{i+1}/{len(phenos_groups)} {pheno_name}")
              #If it's a group of phenos the min function will return the lowest count across all phenos
              tmp_df = df[pheno_list].min(axis =1)
              covs = tmp_df.index[tmp_df >= ~{threshold_cov_count}].tolist()
              missing_covs = [elem for elem in covariates if elem not in covs]
              if missing_covs:tmp_err.write(f"{pheno_name}\t{','.join(missing_covs)}\n")
              o.write(f"{pheno_name}\t{','.join(covs)}\n")
      
      CODE
       >>>
      output {
	File cov_pheno_map = outfile
	File cov_pheno_warning = sub(outfile,'.txt','.err.txt')
      }
  
  runtime {
    cpu: "1"
    docker: "${docker}"
    memory: "${disk_size} GB"
    disks: "local-disk ${disk_size} HDD"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: "1"
  }

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
