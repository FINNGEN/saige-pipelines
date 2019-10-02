task null {
    
    String pheno
    File bedfile
    File bimfile = sub(bedfile, "\\.bed$", ".bim")
    File famfile = sub(bedfile, "\\.bed$", ".fam")
    String prefix = basename(bedfile, ".bed") + "-" + pheno
    File phenofile
    String covariates
    String sampleID
    String traitType
    String docker
    Int cpu
    Int mem
    String loco
    Float traceCVcutoff
    Float ratioCVcutoff
    Int minCovariateCount

    command {

        step1_fitNULLGLMM.R \
            --plinkFile=${sub(bedfile, "\\.bed$", "")} \
            --phenoFile=${phenofile} \
            --phenoCol=${pheno} \
            --covarColList=${covariates} \
            --sampleIDColinphenoFile=${sampleID} \
            --traitType=${traitType} \
            --outputPrefix=${prefix} \
            --nThreads=${cpu} \
            --LOCO=${loco} \
            --traceCVcutoff ${traceCVcutoff} \
            --ratioCVcutoff ${ratioCVcutoff} \
            --minCovariateCount ${minCovariateCount}
    }

    output {

        File modelfile = prefix + ".rda"
        File varianceratio = prefix + ".varianceRatio.txt"
        File randommarkers = prefix + "_30markers.SAIGE.results.txt"
    }

    runtime {

        docker: "${docker}"
        cpu: "${cpu}"
        memory: "${mem} GB"
        disks: "local-disk 20 HDD"
        zones: "europe-west1-b"
        preemptible: 2
        noAddress: true
    }
}

workflow saige_null {

	String docker
    File phenolistfile
    Array[String] phenos = read_lines(phenolistfile)
    String loco

    scatter (pheno in phenos) {

        call null {
            input: docker=docker, pheno=pheno, loco=loco
        }
    }
}
