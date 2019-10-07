task test {

	File nullfile
	File varianceratiofile = sub(nullfile, ".rda", ".varianceRatio.txt")
	File bgenfile
	File samplefile
	Int minmac
	String docker
	Int cpu
    Float mem
    String loco
	String outfile = basename(nullfile, ".rda") + "-" + basename(bgenfile) + ".SAIGE.txt"

    command {

        step2_SPAtests.R \
            --bgenFile=${bgenfile} \
            --minMAC=${minmac} \
            --sampleFile=${samplefile} \
            --GMMATmodelFile=${nullfile} \
            --varianceRatioFile=${varianceratiofile} \
            --SAIGEOutputFile=${outfile} \
            --numLinesOutput=1000 \
            --IsOutputAFinCaseCtrl=TRUE \
            --LOCO=${loco}

    }

	output {
        File out = outfile

    }

    runtime {

        docker: "${docker}"
        cpu: "${cpu}"
        memory: "${mem} GB"
        disks: "local-disk 5 HDD"
        zones: "europe-west1-b"
        preemptible: 2
        noAddress: true
    }
}


task combine {

    String pheno
    Array[File] results
    File blacklist
    Float info_threshold
    String chrcol
    String p_valcol
    String bp_col
    Int loglog_pval
    String docker
    Int cpu
    Float mem

    command <<<

        for file in ${sep=" " results}; do
            if [[ $file == *.gz ]]
            then
                gunzip -c $file > $file"DATAUNZIP"
            else
                mv $file `basename $file`"DATAUNZIP"
            fi
        done

        cat <(head -n 1 `basename ${results[0]}"DATAUNZIP"` | tr ' ' '\t') \
        <(awk 'FNR>1 { printf "%s\t%d\t%s\t%s\t%s\t%s\t%.2f\t%.3e\t%.4f\t%d\t%.4f\t%.4f\t%.4f\t%.3e\t%.3e\t%d\t%.3e\t%.3e\t%.3e\t%.3e\n", \
        $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20 }' \
        `find *DATAUNZIP | sort -V | tr '\n' ' '`) | sort -k1,1V -k2,2g -s > ${pheno} && \
        postscript.py ${pheno} ${blacklist} ${info_threshold} > ${pheno}.pheweb && \
        qqplot.R --file ${pheno}.pheweb --bp_col "${bp_col}" --pval_col "${p_valcol}" --chrcol "${chrcol}" --loglog_pval ${loglog_pval} && \
        bgzip ${pheno} && \
        tabix -S 1 -b 2 -e 2 -s 1 ${pheno}.gz && \
        bgzip ${pheno}.pheweb && \
        tabix -S 1 -b 2 -e 2 -s 1 ${pheno}.pheweb.gz

    >>>

    output {
        File out = pheno + ".gz"
        File out_ind = pheno + ".gz.tbi"
        File out_pheweb = pheno + ".pheweb.gz"
        File out_pheweb_ind = pheno + ".pheweb.gz.tbi"
        File qq = pheno + ".pheweb_" + p_valcol + "_qqplot.png"
        File manh = pheno + ".pheweb_" + p_valcol + "_manhattan.png"
        File manh_loglog = pheno + ".pheweb_" + p_valcol + "_manhattan_loglog.png"
        File quantiles = pheno + ".pheweb_" + p_valcol + "_qquantiles.txt"
    }

    runtime {
        docker: "${docker}"
        cpu: "${cpu}"
        memory: "${mem} GB"
        disks: "local-disk 10 HDD"
        zones: "europe-west1-b"
        preemptible: 2
        noAddress: true
    }
}

workflow test_combine {

    String docker
    String pheno
    String nullfile
    File bgenlistfile
    Array[String] bgenfiles = read_lines(bgenlistfile)
    String loco

    scatter (bgenfile in bgenfiles) {
        call test {
            input: docker=docker, nullfile=nullfile, bgenfile=bgenfile, loco=loco
        }
    }

    call combine {
        input: pheno=pheno, results=test.out, docker=docker
    }
}
