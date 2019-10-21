rule STAR_to_leafcutter_junc:
    input:
        config["RNASeqMapping"]["PathTo_STAR_alignment_Results"] + "{sample}/SJ.out.tab"
    output:
        "sQTL_mapping/juncfiles/{sample}.junc"
    log:
        "logs/sQTL_mapping/STAR_to_leafcutter_junc/{sample}.log"
    shell:
        """
        awk -F'\\t' -v OFS='\\t' '$4==1 && $1!="MT" {{ print $1,$2,$3,".",$7,"+" }} $4==2&& $1!="MT" {{ print $1,$2,$3,".",$7,"-" }}' {input} > {output}
        """

rule make_leafcutter_juncfile:
    input:
        expand ("sQTL_mapping/juncfiles/{sample}.junc", sample=samples["sample"] ),
    output:
        "sQTL_mapping/leafcutter/juncfilelist.txt"
    params:
        SamplesToRemove = config["sQTL_mapping"]["sample_exclude_list"]
    run:
        import os
        SamplesToRemove = open(params.SamplesToRemove, 'r').read().split('\n')
        with open(output[0], "w") as out: 
            for filepath in input:
                samplename = os.path.basename(filepath).split(".junc")[0]
                if samplename not in  SamplesToRemove:
                    out.write(filepath + '\n')

rule leafcutter_cluster:
    input:
        "sQTL_mapping/leafcutter/juncfilelist.txt",
    output:
        "sQTL_mapping/leafcutter/clustering/leafcutter_perind.counts.gz",
        "sQTL_mapping/leafcutter/clustering/leafcutter_perind_numers.counts.gz"
    log:
        "logs/sQTL_mapping/leafcutter_cluster.log"
    shell:
        """
        leafcutter_cluster.py -j {input} -r sQTL_mapping/leafcutter/clustering/ &> {log}
        """

rule leafcutter_prepare_phenotype_table:
    input:
        counts = "sQTL_mapping/leafcutter/clustering/leafcutter_perind.counts.gz",
        blacklist_chromosomes = config["sQTL_mapping"]["chromosome_blacklist"]
    output:
        phenotypes = "sQTL_mapping/leafcutter/clustering/leafcutter_perind.counts.gz.qqnorm_Catted.txt",
        PCs = "sQTL_mapping/leafcutter/clustering/leafcutter_perind.counts.gz.PCs"
    shell:
        """
        ~/miniconda3/bin/python2.7 ~/CurrentProjects/leafcutter/scripts/prepare_phenotype_table.py -p 13 --ChromosomeBlackList {input.blacklist_chromosomes}  {input.counts}

        cat <(head -1 sQTL_mapping/leafcutter/clustering/leafcutter_perind.counts.gz.qqnorm_chr3) <(awk 'FNR>1' sQTL_mapping/leafcutter/clustering/leafcutter_perind.counts.gz.qqnorm_chr*) > {output.phenotypes}
        """

rule prepare_MatrixEQTL_for_sQTL:
    input:
        phenotypes = "sQTL_mapping/leafcutter/clustering/leafcutter_perind.counts.gz.qqnorm_Catted.txt",
        fam = "eQTL_mapping/plink/ForAssociationTesting.fam",
    output:
        phenotypes = "sQTL_mapping/MatrixEQTL/sQTL_phenoTable.txt",
        phenotypesReordered = "sQTL_mapping/MatrixEQTL/sQTL_phenoTable.Reordered.txt",
        intron_locs = "sQTL_mapping/MatrixEQTL/sQTL_intron.locs"
    shell:
        """
        cut -d $'\\t' -f 4- {input.phenotypes} > {output.phenotypes}
        awk -F'\\t' -v OFS='\\t' 'BEGIN {{ print "gene", "chr", "start", "stop" }} NR>1 {{ split($4,a,":"); print $4, a[1], a[2], a[3] }}' {input.phenotypes} > {output.intron_locs}
        Rscript scripts/ReorderPhenotypeTableForMatrixEQTL.R {output.phenotypes} {input.fam} {output.phenotypesReordered}
        """

MatrixSQTLModels = expand("sQTL_mapping/MatrixEQTL/Results/Results.PCs.{B}.covariates.txt", B= list(range (config["sQTL_mapping"]["RNASeqPC_min"], config["sQTL_mapping"]["RNASeqPC_max"] + 1) ))
MatrixSQTLCovariates = expand("sQTL_mapping/Covariates/FromLeafcutter.PCs.{M}.covariates.txt", M=list(range (config["sQTL_mapping"]["RNASeqPC_min"], config["sQTL_mapping"]["RNASeqPC_max"] + 1)))

rule make_covariate_file_sQTL:
    """
    Use R to make covariate file that matches order of fam file
    """
    input:
        fam = "eQTL_mapping/plink/ForAssociationTesting.fam",
        leafcutterPCs = "sQTL_mapping/leafcutter/clustering/leafcutter_perind.counts.gz.PCs",
    output:
        MatrixSQTLCovariates
    params:
        PC_min = config["sQTL_mapping"]["RNASeqPC_min"],
        PC_max = config["sQTL_mapping"]["RNASeqPC_max"],
    shell:
        """
        Rscript scripts/SelectNtoM_PCs_toCovariateFiles.R  {input.leafcutterPCs} {input.fam} sQTL_mapping/Covariates/FromLeafcutter.PCs. {params.PC_min} {params.PC_max}
        """

rule MatrixEQTL_sQTL:
    """Matrix EQTL script performs one cis-eqtl scan with real data and one
    scan with permutated sample labels for phenotypes for an empirical null."""
    input:
        snps = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.snps",
        snp_locs = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.snploc",
        phenotypes = "sQTL_mapping/MatrixEQTL/sQTL_phenoTable.Reordered.txt",
        gene_loc = "sQTL_mapping/MatrixEQTL/sQTL_intron.locs",
        covariates = "sQTL_mapping/Covariates/FromLeafcutter.PCs.{covariate_set}.covariates.txt",
        GRM = "eQTL_mapping/Kinship/GRM.cXX.txt",
    output:
        results = "sQTL_mapping/MatrixEQTL/Results/Results.PCs.{covariate_set}.covariates.txt",
        fig = "sQTL_mapping/MatrixEQTL/Results/images/Results.{covariate_set}.png",
        permutated_fig = "sQTL_mapping/MatrixEQTL/Results/images/PermutatedResults.{covariate_set}.png",
        permuted_results = "sQTL_mapping/MatrixEQTL/Results/PermutatedResults.{covariate_set}.txt",
    log:
        "logs/sQTL_mapping/MatrixEQTL/{covariate_set}.log"
    params:
        cis_window = config["sQTL_mapping"]["cis_window"]
    shell:
        """
        Rscript scripts/MatrixEqtl_Cis.R {input.snps} {input.snp_locs} {input.phenotypes} {input.gene_loc} {input.covariates} {input.GRM} {output.results} {output.fig} {output.permuted_results} {output.permutated_fig} {params.cis_window} &> {log}
        """

rule PlotPCsVsSQTLs:
    input:
        MatrixSQTLModels
    output:
        CattedResult = "sQTL_mapping/MatrixEQTL/Results/ConcatenatedResult.txt",
        Plot = "sQTL_mapping/MatrixEQTL/Results/images/PCsVsSQTLs.pdf"
    log:
        "logs/sQTL_mapping/PlotPCsVsEQTLs.log"
    shell:
        """
        awk -F'\\t' -v OFS='\\t' 'FNR>1 && $6<0.3 {{ print $1,$2,$3,$4,$5,$6,FILENAME  }}' {input} > {output.CattedResult}
        Rscript scripts/Plot_EQTLs_vs_PCs.R {output.CattedResult} {output.Plot}
        """

rule MatrixEQTL_BestModelFromConfigFullResults:
    """
    Matrix eQTL with full output for every snp-gene pair. Also with best p-value for each gene.
    """
    input:
        snps = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.snps",
        snp_locs = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.snploc",
        phenotypes = "sQTL_mapping/MatrixEQTL/sQTL_phenoTable.Reordered.txt",
        gene_loc = "sQTL_mapping/MatrixEQTL/sQTL_intron.locs",
        covariates = "sQTL_mapping/Covariates/FromLeafcutter.PCs.{B}.covariates.txt".format(B=config["sQTL_mapping"]["IdealNumberPC"]),
        GRM = "eQTL_mapping/Kinship/GRM.cXX.txt",
    output:
        results = "sQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/Results.txt",
        BestGenePvals = "sQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/BestPvalueNonPermuted.txt",
        fig = "sQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/images/Results.png",
    log:
        "logs/eQTL_mapping/MatrixEQTL/ConfigCovariateModel.log"
    params:
        cis_window = config["sQTL_mapping"]["cis_window"]
    shell:
        """
        Rscript scripts/MatrixEqtl_Cis.AllPvals.R {input.snps} {input.snp_locs} {input.phenotypes} {input.gene_loc} {input.covariates} {input.GRM} {output.results} {output.fig} {params.cis_window} {output.BestGenePvals} &> {log}
        """

rule splitForEigenMT:
    input:
        snps = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.snps",
        snp_locs = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.snploc",
        results = "sQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/Results.txt",
    output:
        snps = "sQTL_mapping/MatrixEQTL/eigenMT/{chromosome}/gen.txt",
        snp_locs =  "sQTL_mapping/MatrixEQTL/eigenMT/{chromosome}/gen.positions.txt",
        results = "sQTL_mapping/MatrixEQTL/eigenMT/{chromosome}/qtls.txt"
    log:
        "logs/eQTL_mapping/splitForEigenMT/{chromosome}.txt"
    shell:
        """
        awk -F'\\t' -v OFS='\\t' 'NR==1 || /^ID\.{wildcards.chromosome}\./' {input.snps} > {output.snps}
        awk -F'\\t' -v OFS='\\t' 'NR==1 || /^ID\.{wildcards.chromosome}\./' {input.snp_locs} > {output.snp_locs}
        awk -F'\\t' -v OFS='\\t' 'NR==1 || /^ID\.{wildcards.chromosome}\./ {{print $1, $2, $3, $4, $5, $6}}' {input.results} > {output.results}
        """

rule EigenMT:
    # EigenMT is an alternative to permutation testing to get phenotype-level
    # P-values. I used it for identifying eGenes. However I think it was a bit
    # buggy. I had to edit the source code which you can fork from my github.
    # Also, depending on the python executable is pointed to by env (in the
    # shebang line in eigenMT) it might be necessary to install some extra
    # things. Or maybe you can try making a separate conda environment for
    # running eigenMT and using the snakemake conda directive. (I included an
    # environment in this repo but no guarantees that it works out of the box)
    input:
        snps = "sQTL_mapping/MatrixEQTL/eigenMT/{chromosome}/gen.txt",
        snp_locs =  "sQTL_mapping/MatrixEQTL/eigenMT/{chromosome}/gen.positions.txt",
        results = "sQTL_mapping/MatrixEQTL/eigenMT/{chromosome}/qtls.txt",
        gene_loc = "sQTL_mapping/MatrixEQTL/sQTL_intron.locs",
    output:
        "sQTL_mapping/MatrixEQTL/eigenMT/{chromosome}/results.txt"
    log:
        "logs/eQTL_mapping/EigenMT/{chromosome}.log"
    conda:
        "../envs/eigenMT.yaml"
    shell:
        """
        eigenMT.py --QTL {input.results} --GEN {input.snps} --GENPOS {input.snp_locs} --PHEPOS {input.gene_loc} --cis_dist 250000 --OUT {output} --CHROM {wildcards.chromosome} &> {log}
        """

rule catEiginMT:
    input:
        expand("sQTL_mapping/MatrixEQTL/eigenMT/{chromosome}/results.txt", chromosome=set(contigs) - set(blacklist_contigs)),
    output:
        "sQTL_mapping/Chimp.eigenMT.txt.gz",
    log:
        "logs/eQTL_mapping/catEiginMT.log"
    shell:
        """
        (cat <(cat {input} | head -1) <(cat {input} | grep -v '^snps') | gzip - > {output}) &> {log}
        """

rule sQTL_permutations:
    input:
        snps = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.snps",
        snp_locs = "eQTL_mapping/MatrixEQTL/ForAssociationTesting.snploc",
        phenotypes = "sQTL_mapping/MatrixEQTL/sQTL_phenoTable.Reordered.txt",
        gene_loc = "sQTL_mapping/MatrixEQTL/sQTL_intron.locs",
        covariates = "sQTL_mapping/Covariates/FromLeafcutter.PCs.{B}.covariates.txt".format(B=config["sQTL_mapping"]["IdealNumberPC"]),
        GRM = "eQTL_mapping/Kinship/GRM.cXX.txt",
    output:
        results = temp("sQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/Permutations/Chunk.{n}.txt"),
    params:
        InitialSeed = lambda wildcards: int(wildcards.n) * int(config["sQTL_mapping"]["PermutationChunkSize"]),
        NumberPermutations = config["sQTL_mapping"]["PermutationChunkSize"],
        cis_window = config["sQTL_mapping"]["cis_window"]
    log:
        "logs/sQTL_mapping/MatrixEQTL/Permutations/Chunk.{n}.log"
    shell:
        """
        Rscript scripts/MatrixEqtl_Cis_Permutations.R {input.snps} {input.snp_locs} {input.phenotypes} {input.gene_loc} {input.covariates} {input.GRM} {output.results} {params.NumberPermutations} {params.InitialSeed} {params.cis_window} > {log}
        """

LastChunk=int(config["sQTL_mapping"]["NumberPermutationChunks"])-1
rule MergePermutationChunks:
    input:
        Chunks = expand("sQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/Permutations/Chunk.{n}.txt", n=range(0, int(config["sQTL_mapping"]["NumberPermutationChunks"])))
    output:
        "sQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/PermutationsCombined.txt"
    log:
        "sQTL_mapping/MergePermutationChunks.log"
    shell:
        """
        cat sQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/Permutations/Chunk.0.txt > {output}
        for i in {{1..{LastChunk}}}; do
            tail -n +2 sQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/Permutations/Chunk.${{i}}.txt >> {output}
        done
        """

rule PermutationTesting:
    input:
        Permutations_MinPvaluePerTrait = "sQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/PermutationsCombined.txt",
        Actual_MinPvaluePerTrait = "sQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/BestPvalueNonPermuted.txt",
    output:
        PermutationTestResults = "sQTL_mapping/MatrixEQTL/ConfigCovariateModelResults/PermutationTestResults.txt"
    log:
        "logs/sQTL_mapping/PermutationTesting.log"
    shell:
        """
        Rscript scripts/PermutationTest.R {input.Permutations_MinPvaluePerTrait} {input.Actual_MinPvaluePerTrait} {output.PermutationTestResults} &> {log}
        """
