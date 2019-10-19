rule make_plink_file_from_vcf_for_testing:
    input:
        vcf = config["sQTL_mapping"]["vcf"]
    output:
        "eQTL_mapping/plink/Unfiltered.bed"
    log:
        "logs/eQTL_mapping/make_plink_file_from_vcf_for_testing.log"
    shell:
        """
        plink --id-delim '-' --vcf {input.vcf} --vcf-half-call m --allow-extra-chr --make-bed --out eQTL_mapping/plink/Unfiltered &> {log}
        """

rule filter_plink_file_for_testing:
    """
    .fam file output to the gitinclude output directory for convenience, so
    that I can easily access it from my local laptop with git, to edit scripts
    that choose covariates and such with my local RStudio... For the covariate
    files must have samples ordered same as the fam (phenotype) file
    """
    input:
        bed = "eQTL_mapping/plink/Unfiltered.bed"
    output:
        bed = "eQTL_mapping/plink/ForAssociationTesting.bed",
        fam = "eQTL_mapping/plink/ForAssociationTesting.fam",
        bim = "eQTL_mapping/plink/ForAssociationTesting.bim"
    log:
        "logs/eQTL_mapping/filter_plink_file_for_testing.log"
    params:
        keepfam = "--keep-fam " + config["sQTL_mapping"]["keep_fam_list"],
        maf = "--maf " + str(config["sQTL_mapping"]["maf_filter"]),
        remove_ind = "--remove " + config["sQTL_mapping"]["sample_exclude_list"],
        extra = '--geno --hwe 1e-7.5'
    shell:
        """
        plink --bfile eQTL_mapping/plink/Unfiltered  --allow-extra-chr --memory 28000 --make-bed --out eQTL_mapping/plink/ForAssociationTesting {params.keepfam} {params.extra} {params.remove_ind} {params.maf} &> {log}
        """

rule make_plink_file_for_GRM:
    """Note that GRM will ignore individuals with a missing phenotype -9 in the fam file"""
    input:
        bed = "eQTL_mapping/plink/Unfiltered.bed"
    output:
        bed = "eQTL_mapping/Kinship/ForGRM.bed",
        fam = "eQTL_mapping/Kinship/ForGRM.fam"
    params:
        keepfam = "--keep-fam " + config["sQTL_mapping"]["keep_fam_list"],
        maf = "--maf 0.05",
        remove_ind = "--remove " + config["sQTL_mapping"]["sample_exclude_list"],
        extra = '--geno --hwe 1e-7.5'
    shell:
        """
        plink --bfile eQTL_mapping/plink/Unfiltered  --allow-extra-chr --make-bed --out eQTL_mapping/Kinship/ForGRM {params.keepfam} {params.maf} {params.remove_ind} {params.extra}
        sed -i 's/-9$/1/' {output.fam}
        """

rule prune_plink_files_for_GRM:
    input:
        bed = "eQTL_mapping/Kinship/ForGRM.bed",
    output:
        bed =  "eQTL_mapping/Kinship/ForAssociationTesting.pruned.bed",
        fam = "eQTL_mapping/Kinship/ForAssociationTesting.pruned.fam"
    log:
        "logs/eQTL_mapping/prune_plink_files.log"
    shell:
        """
        plink --bfile eQTL_mapping/Kinship/ForGRM --allow-extra-chr --indep-pairwise 50 5 0.5 &> {log}
        plink --bfile eQTL_mapping/Kinship/ForGRM  --allow-extra-chr --extract plink.prune.in --make-bed --out eQTL_mapping/Kinship/ForAssociationTesting.pruned &> {log}
        sed -i 's/-9$/1/' {output.fam}
        rm plink.prune.in plink.prune.out
        """

rule Make_GRM:
    """Genetetic relatedness matrix for gemma LMM"""
    input:
        bed =  "eQTL_mapping/Kinship/ForAssociationTesting.pruned.bed",
    output:
        GRM = "eQTL_mapping/Kinship/GRM.cXX.txt",
    log:
        "logs/eQTL_mapping/make_GRM.log"
    shell:
        """
        gemma -gk 1 -bfile eQTL_mapping/Kinship/ForAssociationTesting.pruned -o GRM -outdir eQTL_mapping/Kinship &> {log}
        """

rule MakeMatrixEQTL_snp_table:
    input:
        plink_bed = "eQTL_mapping/plink/ForAssociationTesting.bed",
    output:
        snps = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.snps",
    params:
        sed_search_delete = config["sQTL_mapping"]["family_name_to_remove"][:-1]
    shell:
        """
        plink --bfile  eQTL_mapping/plink/ForAssociationTesting --allow-extra-chr --recode A-transpose --tab --geno --memory 40000 --out eQTL_mapping/MatrixEQTL/ForAssociationTesting.snps
        perl -lne '/^.+?\\t(.+?\\t).+?\\t.+?\\t.+?\\t.+?\\t(.+$)/ and print "$1$2" ' eQTL_mapping/MatrixEQTL/ForAssociationTesting.snps.traw | sed '1s/_{params}//g' > {output}
        """

rule MakeMatrixEQTL_snp_loc:
    input:
        snps = "eQTL_mapping/plink/ForAssociationTesting.bim",
        plink_bed = "eQTL_mapping/plink/ForAssociationTesting.bed",
    output:
        snp_locs = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.snploc"
    shell:
        """
        awk -F'\\t' -v OFS='\\t' 'BEGIN {{ print "snp","chrom","pos" }} {{ print $2, $1,$4 }}' {input.snps} > {output.snp_locs}
        """

