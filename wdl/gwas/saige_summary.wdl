task summary{
    File input_file
    File finngen_annotation
    File gnomad_annotation
    File finngen_tbi = finngen_annotation + ".tbi"
    File gnomad_tbi = gnomad_annotation + ".tbi"

    Float pval_thresh
    String pheno
    String output_name=pheno+"_summary"
    String docker

    command <<<
        set -euxo pipefail
        python3 <<EOF
        #read file
        fname = "${input_file}"
        finngen_annotation_file = "${finngen_annotation}"
        gnomad_annotation_file  = "${gnomad_annotation}"
        sig_threshold = ${pval_thresh}

        output_name = "${output_name}"

        import pysam
        import gzip
        from typing import NamedTuple

        FGAnnotation = NamedTuple('FGAnnotation',[('gene',str),('consequence',str)])
        GnomadAnnotation = NamedTuple('GnomadAnnotation',[('finaf',str),('nfeaf',str),('rsid',str)])

        def get_header(reader,file):
            with reader(file, "rt") as f:
                l = f.readline()
                return l.strip("\n").split("\t")

        def get_fg_annotation(iterator,variant,gene_idx, consequence_idx,variant_idx) -> FGAnnotation:
            for v in iterator:
                cols = v.strip("\n").split("\t")
                if cols[variant_idx] == variant:
                    return FGAnnotation(cols[gene_idx],cols[consequence_idx])
            return FGAnnotation("","")

        def get_gnomad_annotation(iterator,cpra,finaf_idx, nfeaf_idx, rsid_id,gnomad_cpra) -> GnomadAnnotation:
            for v in iterator:
                cols = v.strip("\n").split("\t")

                if (cols[gd_cpra[0]] == "chr"+cpra[0]) and (cols[gd_cpra[1]] == str(cpra[1])) and (cols[gd_cpra[2]] == cpra[2]) and (cols[gd_cpra[3]] == cpra[3]):
                    return GnomadAnnotation(cols[finaf_idx],cols[nfeaf_idx],cols[rsid_idx])
            return GnomadAnnotation(".",".",".")

        #required columns
        fg_req_cols=["#variant","gene_most_severe","most_severe"]
        gd_req_cols=["fin.AF","nfsee.AF","enrichment_nfsee","rsid"]
        #open finngen annotation tabix
        fg_tabix = pysam.TabixFile(finngen_annotation_file,parser=None)
        #get fg header column positions
        fg_header = get_header(gzip.open, finngen_annotation_file)
        fg_idx = {v:i for (i,v) in enumerate(fg_header)}
        #open gnomad annotation tabix
        gnomad_tabix = pysam.TabixFile(gnomad_annotation_file,parser=None)
        gd_header = get_header(gzip.open, gnomad_annotation_file)
        gd_idx = {v:i for (i,v) in enumerate(gd_header)}
        gd_cpra = (
            gd_idx["chrom"],
            gd_idx["pos"],
            gd_idx["ref"],
            gd_idx["alt"]
        )
        #check for fg column existence
        if not all([a in fg_header for a in fg_req_cols]):
            raise Exception("Not all columns present in FinnGen annotation! Aborting...")
        #check for gnomad column existence
        if not all([a in gd_header for a in gd_req_cols]):
            raise Exception("Not all columns present in Gnomad annotation! Aborting...")
        var_idx, gene_idx, cons_idx = (fg_idx["#variant"],fg_idx["gene_most_severe"],fg_idx["most_severe"])

        finaf_idx, nfeaf_idx, rsid_idx = (gd_idx["fin.AF"],gd_idx["nfsee.AF"],gd_idx["rsid"])



        with gzip.open(fname, "rt") as file:
            #open output file
            with open(output_name,"w") as outfile:
                #read header
                header = file.readline().strip("\n").split('\t')
                #find p-value index
                header_idx = {v:i for (i,v) in enumerate(header)}
                pval_idx = header_idx["pval"]
                cid,pid,rid,aid = (
                    header_idx["#chrom"],
                    header_idx["pos"],
                    header_idx["ref"],
                    header_idx["alt"]
                )

                #add gene name, consequence, gnomad finnish af, nfsee af, rsid
                header.extend(["fin.AF","nfsee.AF","rsid","gene_most_severe","most_severe"])
                outfile.write("\t".join(header)+"\n")
                #read lines
                for line in file:
                    line_columns = line.strip("\n").split('\t')
                    pvalue = float(line_columns[pval_idx])
                    if pvalue < sig_threshold:
                        cpra= (line_columns[cid],int(float(line_columns[pid])),line_columns[rid],line_columns[aid])
                        variant = "{}:{}:{}:{}".format(cpra[0],cpra[1],cpra[2],cpra[3])
                        fg_c = cpra[0].replace("chr","").replace("X","23").replace("Y","24").replace("M","25").replace("MT","25")
                        gd_c = "chr" + cpra[0].replace("chr","").replace("23","X").replace("24","Y").replace("25","M")
                        #annotate
                        fg_iter = fg_tabix.fetch(fg_c,cpra[1]-1, cpra[1])
                        fg_a = get_fg_annotation(fg_iter,variant, gene_idx,cons_idx,var_idx)

                        #annotate
                        gnomad_iter = gnomad_tabix.fetch(gd_c,cpra[1]-1, cpra[1])
                        gd_a = get_gnomad_annotation(gnomad_iter,cpra,finaf_idx, nfeaf_idx, rsid_idx,gd_cpra)

                        line_columns.extend([
                            gd_a.finaf,
                            gd_a.nfeaf,
                            gd_a.rsid,
                            fg_a.gene,
                            fg_a.consequence,
                        ])
                        #gather row
                        #write to file
                        outfile.write("\t".join(line_columns)+"\n")
        print("summary created")
        EOF
        bgzip ${output_name}

    >>>

    output {
        File out = output_name+".gz"
    }

    runtime {
        docker: "${docker}"
        cpu: 1
        memory: "10 GB"
        disks: "local-disk 200 HDD"
        zones: "europe-west1-b"
        preemptible: 2
    }
}
