
######################################################
"""
            This rule is used for cleaning the data extracted from the TCGA database
"""
#######################################################

##### unzip all the data and merge and rename the expression data to two files: FPKM and COUNTS

rule TCGA_unzip_and_merge_expression_data:
    input:
        TCGA_dir = "data/TCGA"

    output:
        fpkm = "data/TCGA/expression_FPKM.tsv",
        raw_count = "data/TCGA/expression_count.tsv",

    shell:
        """
        mkdir {input.TCGA_dir}/expression &&
        TAR_FILES=`find {input.TCGA_dir} -name *.gz | grep "gdc_download"` &&
        for TAR_FILE in $TAR_FILES ; do
            tar -C  {input.TCGA_dir}/expression -xf $TAR_FILE
        done &&
        rm  {input.TCGA_dir}/expression/MANIFEST.txt &&

        GZ_FILES=`find {input.TCGA_dir}/expression -name *.gz` &&
        for GZ_FILE in $GZ_FILES ; do
            gunzip $GZ_FILE
        done &&

        mkdir {input.TCGA_dir}/expression2 &&
        COUNTS_FILES=`find {input.TCGA_dir}/expression -name *.counts`
        for COUNTS_FILE in $COUNTS_FILES ; do
            mv $COUNTS_FILE  {input.TCGA_dir}/expression2
        done &&
        FPKM_FILES=`find {input.TCGA_dir}/expression -name *.FPKM.txt`
        for FPKM_FILE in $FPKM_FILES ; do
            mv $FPKM_FILE {input.TCGA_dir}/expression2
        done &&

        FILE_NAMES=`ls {input.TCGA_dir}/expression2`
        for FILE_NAME in $FILE_NAMES ; do
            echo $FILE_NAME
            FILE_PATH={input.TCGA_dir}/expression2/$FILE_NAME

            echo "genes\t$FILE_NAME" > {input.TCGA_dir}/temporary_file

            (cat {input.TCGA_dir}/temporary_file && cat $FILE_PATH) > {input.TCGA_dir}/expression2/test &&
            mv {input.TCGA_dir}/expression2/test $FILE_PATH
        done &&

        paste {input.TCGA_dir}/expression2/*.FPKM.txt > {output.fpkm} &&
        paste {input.TCGA_dir}/expression2/*.counts > {output.raw_count} &&

        rm -r {input.TCGA_dir}/expression &&
        rm -r {input.TCGA_dir}/expression2 &&
        rm {input.TCGA_dir}/temporary_file

        """

######################################################

#### Rule for extract and merge clinical data and id data

rule TCGA_extract_and_merge_clinical_data:
    input:
        TCGA_dir = "data/TCGA",

    output:
        TCGA_id_samples = "data/TCGA/id_samples.tsv",
        TCGA_clinical_data = "data/TCGA/TCGA_clinical_data.tsv"

    shell:
        """
        LGG_file=`find {input.TCGA_dir} -name *gdc_sample_sheet.2020-01-22.tsv | grep "LGG"`
        GBM_file=`find {input.TCGA_dir} -name *gdc_sample_sheet.2020-01-22.tsv | grep "GBM"`
        (cat $LGG_file && tail -n +2 $GBM_file) > {output.TCGA_id_samples} &&

        clinical_TAR_file_LGG=`find {input.TCGA_dir} -name *clinical.cart.2020-01-22.tar.gz | grep "LGG"` &&
        clinical_TAR_file_GBM=`find {input.TCGA_dir} -name *clinical.cart.2020-01-22.tar.gz | grep "GBM"` &&
        tar -C {input.TCGA_dir}/LGG -xf $clinical_TAR_file_LGG &&
        tar -C {input.TCGA_dir}/GBM -xf $clinical_TAR_file_GBM

        clinical_file_LGG=`find {input.TCGA_dir} -name *clinical.tsv | grep "LGG"` &&
        clinical_file_GBM=`find {input.TCGA_dir} -name *clinical.tsv | grep "GBM"` &&
        (cat $clinical_file_LGG && tail -n +2 $clinical_file_GBM) > {output.TCGA_clinical_data}
        """




#######################################################

##### Rule cleaning the expression data and the clinical data

rule TCGA_cleaning:
    input:
        expression_data_TCGA_count = "data/TCGA/expression_count.tsv",
        expression_data_TCGA_fpkm = "data/TCGA/expression_FPKM.tsv",
        id_data_TCGA = "data/TCGA/id_samples.tsv",
        clinical_data_TCGA = "data/TCGA/TCGA_clinical_data.tsv",
        # id_data_TCGA = rules.TCGA_extract_and_merge_clinical_data.output.TCGA_id_samples,
        # clinical_data_TCGA = rules.TCGA_extract_and_merge_clinical_data.output.TCGA_clinical_data,
        # expression_data_TCGA_fpkm = rules.TCGA_unzip_and_merge_expression_data.output.fpkm,
        # expression_data_TCGA_count = rules.TCGA_unzip_and_merge_expression_data.output.raw_count,





    output:
        TCGA_clinical_data_cancer = "data/clinical_data/TCGA_clinical_data_cancer",
        TCGA_clinical_data_normal = "data/clinical_data/TCGA_clinical_data_normal",
        TCGA_expression_data_fpkm_cancer = "data/expression_data/TCGA_expression_data_fpkm_cancer",
        TCGA_expression_data_count_cancer = "data/expression_data/TCGA_expression_data_count_cancer",
        TCGA_expression_data_fpkm_normal = "data/expression_data/TCGA_expression_data_fpkm_normal",
        TCGA_expression_data_count_normal = "data/expression_data/TCGA_expression_data_count_normal",

    script:
        "../../scripts/python_modules/TCGA_cleaning.py"

#######################################################
