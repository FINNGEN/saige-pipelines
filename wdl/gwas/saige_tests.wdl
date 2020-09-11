import "saige_sub.wdl" as sub

workflow saige_tests {

    String docker
    File phenolistfile
    String traitType
    Array[String] phenos = read_lines(phenolistfile)
    String nullfiletemplate
    String loco
    String analysisType

    scatter (pheno in phenos) {

        String nullfile = sub(nullfiletemplate, "PHENO", pheno)

        call sub.test_combine {
            input: docker=docker, pheno=pheno, traitType=traitType, nullfile=nullfile, loco=loco, analysisType=analysisType
        }
    }
}
