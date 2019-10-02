import "saige_sub.wdl" as sub

workflow saige_tests {

	String docker
    File phenolistfile
    Array[String] phenos = read_lines(phenolistfile)
    String nullfiletemplate
    String loco

    scatter (pheno in phenos) {

        String nullfile = sub(nullfiletemplate, "PHENO", pheno)

        call sub.test_combine {
            input: docker=docker, pheno=pheno, nullfile=nullfile, loco=loco
        }
    }
}
