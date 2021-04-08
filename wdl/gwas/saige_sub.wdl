task test {

    File nullfile
    File varianceratiofile = sub(nullfile, ".rda", ".varianceRatio.txt")
    String outfileprefix = basename(nullfile, ".rda") + "-"
    Array[File] bgenfiles
    File samplefile
    Int minmac
    String docker
    String loco
    String analysisType
    Boolean logP
    String logPStr = if logP then "TRUE" else "FALSE"

    command {
        set -euxo pipefail
        python3 <<EOF
        import os
        import subprocess
        import time
        processes = set()
        cmd_prefix = 'export MKL_NUM_THREADS=1; export MKL_DYNAMIC=false; export OMP_NUM_THREADS=1; export OMP_DYNAMIC=false; \
            step2_SPAtests.R \
            --minMAC=${minmac} \
            --sampleFile=${samplefile} \
            --GMMATmodelFile=${nullfile} \
            --varianceRatioFile=${varianceratiofile} \
            --numLinesOutput=1000 \
            --IsOutputAFinCaseCtrl=TRUE \
            --LOCO=${loco} \
            --IsOutputlogPforSingle=${logP} \
            --analysisType=${analysisType} '
        for file in '${sep=" " bgenfiles}'.split(' '):
            cmd = cmd_prefix + '--bgenFile=' + file
            cmd = cmd + ' --SAIGEOutputFile=${outfileprefix}' + os.path.basename(file) + '.SAIGE.txt'
            logfile = open('SAIGE_log_${outfileprefix}' + os.path.basename(file) + '.txt', 'w')
            processes.add(subprocess.Popen(cmd, shell=True, stdout=logfile))
        print(time.strftime("%Y/%m/%d %H:%M:%S") + ' ' + str(len(processes)) + ' processes started', flush=True)
        n_rc0 = 0
        while n_rc0 < len(processes):
            time.sleep(60)
            n_rc0 = 0
            for p in processes:
                p_poll = p.poll()
                if p_poll is not None and p_poll > 0:
                    raise Exception('subprocess returned ' + str(p_poll))
                if p_poll == 0:
                    n_rc0 = n_rc0 + 1
            print(time.strftime("%Y/%m/%d %H:%M:%S") + ' ' + str(n_rc0) + ' processes finished', flush=True)
        EOF
    }

	output {
        Array[File] out = glob("*.SAIGE.txt")
        Array[File] logs = glob("SAIGE_log_*.txt")
    }

    runtime {
        docker: "${docker}"
        cpu: length(bgenfiles)
        memory: (4 * length(bgenfiles)) + " GB"
        disks: "local-disk " + (length(bgenfiles) * ceil(size(bgenfiles[0], "G")) + 5) + " HDD"
        zones: "europe-west1-b"
        preemptible: 2
        noAddress: true
    }
}

task combine {

    String pheno
    String traitType
    Array[Array[File]] results2D
    Array[File] results = flatten(results2D)
    String chrcol
    String p_valcol
    String bp_col
    Int loglog_pval
    String docker
    String prefix
    Boolean logP
    String logPStr = if logP then "True" else "False"

    command <<<
        set -euxo pipefail

        echo "`date` concatenating results to ${prefix}${pheno}.saige.gz"
        cat \
        <(head -n 1 ${results[0]}) \
        <(for file in ${sep=" " results}; do tail -n+2 $file; done) \
        |awk 'NR == 1; NR > 1 {print $0 | "sort -k 1,1V -k 2,2g"}' \
        |awk 'NR==1 { for(i=1;i<=NF;i++) { h[$i]=i }; if (! "POS" in h) { print "ERROR: POS not in saige file header"; exit 1};print } NR>1{ $h["POS"]=sprintf("%d",$h["POS"]);print $0 }' \
        |tr ' ' '\t' \
        |bgzip > ${prefix}${pheno}.saige.gz
        tabix -s1 -b2 -e2 -S1 ${prefix}${pheno}.saige.gz

        echo "`date` converting results to ${prefix}${pheno}.gz"
        python3 <<EOF | sort -k 1,1g -k 2,2g | bgzip > ${prefix}${pheno}.gz
        import math, gzip
        from collections import OrderedDict
        from functools import reduce

        def conv_base(logp):
            return logp / math.log(10)

        def red(obj, func_list):
            return "NA" if obj == "NA" else reduce(lambda o, func: func[0](o, *func[1]) if func[0] is not str.format else func[0](func[1], o), func_list, obj)

        mapping = OrderedDict([
            ("#chrom", ("CHR", [(str.replace, ("chr", "")), (str.replace, ("X", "23")), (str.replace, ("Y", "24")), (str.replace, ("MT", "25")), (str.replace, ("M", "25"))])),
            ("pos", ("POS", [(float, ()), (int, ())])), # convert 1e6 to 1000000
            ("ref", ("Allele1", [])),
            ("alt", ("Allele2", [])),
            ("pval", ("p.value", [(float, ()), (math.exp, ()), (str.format, ("{:.2e}"))])) if ${logPStr} else ("pval", ("p.value", [(float, ()), (str.format, ("{:.2e}"))])),
            ("mlogp", ("p.value", [(float, ()), (abs, ()), (conv_base, ()), (str.format, ("{:.4f}"))])) if ${logPStr} else ("mlogp", ("p.value", [(float, ()), (math.log10, ()), (abs, ()), (str.format, ("{:.4f}"))])),
            ("beta", ("BETA", [(float, ()), (str.format, ("{:.5f}"))])),
            ("sebeta", ("SE", [(float, ()), (str.format, ("{:.5f}"))])),
            ("af_alt", ("AF_Allele2", [(float, ()), (str.format, ("{:.2e}"))])),
            ("af_alt_cases", ("AF.Cases", [(float, ()), (str.format, ("{:.2e}"))])),
            ("af_alt_controls", ("AF.Controls", [(float, ()), (str.format, ("{:.2e}"))])),
            ("n_hom_cases", ("homN_Allele2_cases", [(float, ()), (str.format, ("{:.2f}"))])),
            ("n_het_cases", ("hetN_Allele2_cases", [(float, ()), (str.format, ("{:.2f}"))])),
            ("n_hom_ref_cases", ("ref_homN_Allele2_cases", [(float, ()), (str.format, ("{:.2f}"))])),
            ("n_hom_controls", ("homN_Allele2_ctrls", [(float, ()), (str.format, ("{:.2f}"))])),
            ("n_het_controls", ("hetN_Allele2_ctrls", [(float, ()), (str.format, ("{:.2f}"))])),
            ("n_hom_ref_controls", ("ref_homN_Allele2_ctrls", [(float, ()), (str.format, ("{:.2f}"))])),
        ])

        with gzip.open("${prefix}${pheno}.saige.gz", 'rt') as f:
            header = {h:i for i,h in enumerate(f.readline().strip().split('\t'))}
            for col in [v[0] for v in mapping.values()]:
                if col not in header.keys():
                    raise Exception('column ' + col + ' not in given file')
            print('\t'.join(mapping.keys()))
            for line in f:
                s = line.strip().split('\t')
                print('\t'.join(str(red(s[header[v[0]]], v[1])) for v in mapping.values()))
        EOF
        echo "`date` plotting qq and manha"
        qqplot.R --file ${prefix}${pheno}.gz --bp_col "${bp_col}" --pval_col "${p_valcol}" --chrcol "${chrcol}" --loglog_pval ${loglog_pval}
        echo "`date` tabixing"
        tabix -S 1 -b 2 -e 2 -s 1 ${prefix}${pheno}.gz
        echo "`date` done"
    >>>

    output {
        File saige_out = prefix + pheno + ".saige.gz"
        File saige_out_ind = prefix + pheno + ".saige.gz.tbi"
        File out = prefix + pheno + ".gz"
        File out_ind = prefix + pheno + ".gz.tbi"
        Array[File] pngs = glob("*.png")
        Array[File] quantiles = glob("*qquantiles.txt")
    }

    runtime {
        docker: "${docker}"
        cpu: 1
        memory: "10 GB"
        disks: "local-disk 200 HDD"
        zones: "europe-west1-b"
        preemptible: 2
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
    String analysisType

    scatter (bgenfiles in bgenfiles2D) {
        call test {
            input: docker=docker, nullfile=nullfile, bgenfiles=bgenfiles, loco=loco, analysisType=analysisType
        }
    }

    call combine {
        input: pheno=pheno, traitType=traitType, results2D=test.out, docker=docker
    }

    output {
        File saige_out = combine.saige_out
        File saige_out_ind = combine.saige_out_ind
        File out = combine.out
        File out_ind = combine.out_ind
        Array[File] pngs = combine.pngs
        Array[File] quantiles = combine.quantiles
    }
}
