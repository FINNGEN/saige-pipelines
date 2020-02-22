task test {

	File nullfile
	File varianceratiofile = sub(nullfile, ".rda", ".varianceRatio.txt")
    String outfileprefix = basename(nullfile, ".rda") + "-"
	Array[File] bgenfiles
	File samplefile
	Int minmac
	String docker
    String loco

    command {

        python3 <<EOF
        import os
        import shlex
        import subprocess
        processes = set()
        cmd_prefix = 'step2_SPAtests.R \
            --minMAC=${minmac} \
            --sampleFile=${samplefile} \
            --GMMATmodelFile=${nullfile} \
            --varianceRatioFile=${varianceratiofile} \
            --numLinesOutput=1000 \
            --IsOutputAFinCaseCtrl=TRUE \
            --LOCO=${loco} '
        for file in '${sep=" " bgenfiles}'.split(' '):
            cmd = cmd_prefix + '--bgenFile=' + file
            cmd = cmd + ' --SAIGEOutputFile=${outfileprefix}' + os.path.basename(file) + '.SAIGE.txt'
            logfile = open('SAIGE_log_${outfileprefix}' + os.path.basename(file) + '.txt', 'w')
            processes.add(subprocess.Popen(shlex.split(cmd), stdout=logfile))
        for p in processes:
            if p.poll() is None:
                p.wait()
        EOF
    }

	output {
        Array[File] out = glob("*.SAIGE.txt")
        Array[File] logs = glob("SAIGE_log_*.txt")
    }

    runtime {
        docker: "${docker}"
        cpu: length(bgenfiles)
        memory: (3.75 * length(bgenfiles)) + " GB"
        disks: "local-disk " + (length(bgenfiles) * ceil(size(bgenfiles[0], "G")) + 1) + " HDD"
        zones: "europe-west1-b"
        preemptible: 0
        noAddress: true
    }
}


task combine {

    String pheno
    String traitType
    Array[Array[File]] results2D
    Array[File] results = flatten(results2D)
    File blacklist
    Float info_threshold
    String chrcol
    String p_valcol
    String bp_col
    Int loglog_pval
    String docker

    command <<<

        for file in ${sep=" " results}; do
            if [[ $file == *.gz ]]
            then
                gunzip -c $file > $file"DATAUNZIP"
            else
                mv $file `basename $file`"DATAUNZIP"
            fi
        done

        if [[ ${traitType} == "binary" ]]; then
            cat <(head -n 1 `basename ${results[0]}"DATAUNZIP"` | tr ' ' '\t') \
            <(awk 'FNR>1 { printf "%s\t%d\t%s\t%s\t%s\t%s\t%.2f\t%.3e\t%.2f\t%.2f\t%.4f\t%d\t%.4f\t%.4f\t%.4f\t%.3e\t%.3e\t%d\t%.3e\t%.3e\t%.3e\t%.3e\n", \
            $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22 }' \
            `find *DATAUNZIP | sort -V | tr '\n' ' '`) | sort -k1,1V -k2,2g -s > ${pheno}
        else
            cat <(head -n 1 `basename ${results[0]}"DATAUNZIP"` | tr ' ' '\t') \
            <(awk 'FNR>1 { printf "%s\t%d\t%s\t%s\t%s\t%s\t%.2f\t%.3e\t%.2f\t%.2f\t%.4f\t%d\t%.4f\t%.4f\t%.4f\t%.3e\t%.3e\t%d\t%.3e\t%.3e\n", \
            $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20 }' \
            `find *DATAUNZIP | sort -V | tr '\n' ' '`) | sort -k1,1V -k2,2g -s > ${pheno}
        fi
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
        cpu: 1
        memory: "7 GB"
        disks: "local-disk 20 HDD"
        zones: "europe-west1-b"
        preemptible: 0
        noAddress: true
    }
}

workflow test_combine {

    String docker
    String pheno
    String traitType
    String nullfile
    File bgenlistfile
    Array[Array[String]] bgenfiles2D = read_tsv(bgenlistfile)
    String loco

    scatter (bgenfiles in bgenfiles2D) {
        call test {
            input: docker=docker, nullfile=nullfile, bgenfiles=bgenfiles, loco=loco
        }
    }

    call combine {
        input: pheno=pheno, traitType=traitType, results2D=test.out, docker=docker
    }
}
