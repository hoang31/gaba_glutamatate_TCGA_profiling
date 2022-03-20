


# **GABA and Glutamate Transcriptomic Profiling in glioma**

## **1st step : install miniconda**
You will see all the documentation for installing in the main website (https://docs.conda.io/en/latest/miniconda.html).

## **2nd step : create snakemake environment**

    conda create --name smake_chapter2 -c bioconda -c conda-forge snakemake=5.10.0
    conda activate smake_chapter2
    conda install -c anaconda pandas=0.25.3 numpy


## **3th step : download all data**

In this pipeline, some data will be downloaded automatically. These data are :
  - CGGA1 (expression and clinical data)
  - CGGA2 (expression and clinical data)

For the other database :
  - TCGA : we need to download all the data in the TCGA portal (https://portal.gdc.cancer.gov/). For each project (LGG and GBM projects), we need the biospecimen, clinical, sample sheet, meta data and expression data available in the cart.



The data strucure to follow for the snakemake pipeline :  
---> PROJECT directory  
------> snakefile  
------> data  
----------> TCGA  
----------> CGGA1  
----------> CGGA2  

## **4th step : Execute the snakemake pipeline**
    snakemake --use-conda
