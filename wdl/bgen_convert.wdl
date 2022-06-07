version 1.0

workflow convert_bgen {


  input{
    Array[String] chrom_list 
    String docker
    String name
  }
  scatter (chrom in chrom_list){
    
    call chrom_convert {
      input :
      chrom = chrom,
      docker = docker,
      name = name
      }     
     }
          
}

task chrom_convert {
  input {
    String chrom
    String name
    String docker
    String chromPath
    String bargs
    Int disk_factor
    File variants

  }
  # get path to vcf and tabix file
  
  File cFile = sub(chromPath,"CHROM",chrom)
  String name_chrom =  name + "_" + chrom
  File tbiFile = cFile + '.tbi'
  # get chrom size and convert to disk size & runtime options
  Int chrom_size = ceil(size(cFile,"GB"))
  String disk_size =disk_factor * chrom_size +2
  #bgen conversion options 
  
  # CPU AND MEM
  Int tmp_cpu =  if chrom_size < 48 then chrom_size*2 else 95
  Int cpu = if tmp_cpu < 8 then 8 else tmp_cpu #cpus in case of splitting
  
  Int tmp_mem = if chrom_size < 64 then chrom_size else 64
  Int mem = if tmp_mem < 16 then 16 else tmp_mem
   
  Int chrom_int = chrom
  
  command <<<
    cat ~{variants} | grep -w ~{chrom} | cut -f2 > ./variants.txt
    echo ~{disk_size}  ~{mem} ~{cpu} && df -h
    
    python3 /Scripts/annotate.py \
    --cFile ~{cFile} \
    --tbiFile ~{tbiFile} \
    --oPath "/cromwell_root/" \
    --vcf-variants ./variants.txt \
    --name ~{name_chrom} \
    --annotate \
    --split \
    -b  \
    --bargs ~{bargs} \
    | tee  chrom_convert_~{name_chrom}.log        
    df -h >> chrom_convert_~{name_chrom}.log

    >>>

    runtime {
        docker: "${docker}"
        cpu: "${cpu}"
        disks: "local-disk " + "${disk_size}" + " HDD"
        bootDiskSizeGb: 20
        memory: "${mem}" + " GB"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
	preemptible: 0
    }

    output {    
        File out_bgen = "/cromwell_root/${name_chrom}/${name_chrom}.bgen"
        File out_bgen_sample = "/cromwell_root/${name_chrom}/${name_chrom}.bgen.sample"
        File out_bgen_index = "/cromwell_root/${name_chrom}/${name_chrom}.bgen.bgi"
        File out_chrom_convert_log  = "chrom_convert_${name_chrom}.log"
	Array[File] logs = glob("/cromwell_root/${name_chrom}/logs/*")
    }
}

