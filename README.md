Introduction


# 2 parts are included in this project
## 1)Models generation
Model generation is done manually on jupyter notebook and is not written into pipeline
## 2)Fragment ranking and protein quantify algrithm

## Pipelines
Useful workflow includes the generalized OpenSWATH workflow and the optimized OpenSWATH workflow which combine OpenSWATH, pyProphet and a python package which used to combine classifier and q value calculation algorithm(gscore) to score fragment ions.

## Main modules
### The generalized OpenSWATH workflow
OpenSWATH + pyProphet + TRIC

### the optimized OpenSWATH workflow
OpenSWATH + pyProphet
Generated classifiers
python package that used to quantify protein at fragment level


The workflows were written using snakemake language with some python language inside, it can be used to proteomics research, especially for DIA MS collected data. It allows users to generate library and get EncyclopeDIA quantitative results.

In order to execute snakemake workflow, the installation of java,python,docker and snakemake is required. The absolute path of the files to generate the library and to get quantitative results should be pass to the workflow as variable values via the command line.

Environment and software supports
openswath
pyprophet
python
snakemake
