# Introduction

This project further optimizes the openSWATH workflow with FDR control at the fragment ion level using machine learning, and quantifies DIA data using different precursor quantification methods.
The provided code contains 1 generalized openSWATH workflow, an optimized DIA data analysis workflow, a classifier training workflow, and a python package for FDR control at the fragment ion level.
The optimized DIA data analysis consists of 3 steps, 1. The discriminant score is calculated from the 4 classifiers trained, and the discriminant score is used to calculate the q value, and then the fragment ions are controlled at the 1% FDR level. 2. The remaining fragment ions are used to calculate the intensities of the precursor ions, mean, topN, median or sum are optional methods. 3. Each run is processed individually and then combined into a quantitative matrix.  

## Files and data

```bash
.
├── data
│   ├── 220126_oswag_library.optimized.decoy.pqp    # Spectra library
│   ├── 220126_oswag_rt_peptides.optimized.pqp	# iRT library
│   ├── bruduer_HFX_44win.txt   # swath window file
│   └── openms.sif  # OpenSWATH singularity container
├── dia_quantification
│   ├── classifier_filter_fragments # Portable python package that can quantify
│   │   ├── filter.py   pyProphet results at fragment ion level
│   │   ├── global_quantify.py
│   │   ├── merge_osw.py
│   │   ├── metrics.py
│   │   ├── predict.py
│   │   └── preprocessing.py
│   └── snakemake
│       ├── classifier_quant.smk    # Optimized pipeline 
│       ├── classifier_quant.yml    # Config file for optimized pipeline
│       ├── pyprophet_feed.smk  # Standard pipeline for baseline results
│       ├── pyprophet_model.smk # Standard pipeline for model generation
│       └── pyprophet.yml   # Config file for standard pipeline
└── model_generation    # Machine learning model and training process
    ├── models  # Trained models
    │   ├── initial.pkl
    │   ├── lnl.pkl
    │   ├── precision_90.pkl
    │   ├── qvalue_cutoff.pkl
    │   ├── standardscaler_90.pkl
    │   ├── standardscaler.pkl  # For initial and lnl classifiers
    │   └── standardscaler_qvalue.pkl
    ├── build_classifiers.ipynb # Model training process
    └── training_data   # Training data preprocessing
        ├── classifier_training_data.smk
        ├── classifier_training_data.yml
        ├── merge_csv_files.py
        └── oswSQL2csv.py
```



## 1. Generalized workflow 

It is the routine workflow of  analyzing DIA data using OpenSwath.  It contains 3 modules:

**① OpenSwath** is used for targeted extraction of chromatogram; 

**②pyProphet** is used for statistical scoring; 

**③TRIC** is used for alignment of multiple runs to generate quantitative matrix. 

It contains 2 pipelines written with Snakemake:

**① pyprophet_model.smk** conducts openSWATH and part of pyprophet workflow. The input is mzml files, swath window file, spectra library and iRT library, the output is pyprophet score model "score_model.osw". This step is also a must for optimized workflow.

**② pyprophet_feed.smk** conducts part of pyprophet workflow, it takes osw files generated from openSWATH and trained scoring_model.osw as input, and generate quantitative matrix.

***pyprophet.yml** is the configure file of generalized workflow.

<img src="https://github.com/Lina0125/QuantifyAtFragmentLevelWithClassifiers/blob/main/data/imgs/generalized.png" alt="generalized" title="The stantard workflow" style="zoom30%;" >

### Environment and software supports

OpenSWATH(singularity container)

pyProphet

TRIC

python

snakemake

### Inputs

*.mzml

swath window file

spectra library

iRT library

### Usage

```bash
$snakemake -s pyprophet_model.smk --configfile pyprophet.yml --core 50
$snakemake -s pyprophet_feed.smk --configfile pyprophet.yml --core 50
```



## 2) Classifiers generation

To keep track of where classifiers came from, a training script used to generate results in the report is provided(build_classifiers.ipynb). The training data for building the classifier comes from the output of OpenSWATH. OpenSWATH extracts the features of transition groups from the chromatogram and stores them in the osw tables. The 12 sub-scores are used as classifying features, and the label is the target and the decoy added to the search space. The reference classifier is trained using SGD first to provide reference for the improvements. The model then improved by precision/recall trade-off, FDR control and confident learning. This script finally trained 4 models and it's corresponding standard scalers.

<img src="https://github.com/Lina0125/QuantifyAtFragmentLevelWithClassifiers/blob/main/data/imgs/model_generation.png" alt="model_generation" style="zoom:30%;" />

 

## 3) Optimized workflow

classifier_quant.smk is written in snakemake and it contains 3 modules: 

1.OpenSWATH for chromatographic extraction

2.pyProphet for statistical analysis data on protein and peptides level

3.classifier_filter_fragments python package first scores every fragment ion for precursors using generated classifiers, then the discriminate score is used to calculate q values for FDR control using gscore python package, after control FDR on fragment ion level, the optimized fragment ion population is used to quantify peptides intensity by mean, median, topN and sum. Every run is processed separately then are formed into a quantitative matrix.

<img src="https://github.com/Lina0125/QuantifyAtFragmentLevelWithClassifiers/blob/main/data/imgs/application.png" alt="application" style="zoom:30%; margin:auto; float:center" />



### Environment and software supports

OpenSWATH(singularity container)

pyProphet

python

python package: classifier_filter_fragments; gscore; numba

### Inputs

*.mzml

swath window file

spectra library

iRT library

classifiers

### Usage:

```bash
$snakemake -s pyprophet_model.smk --configfile pyprophet.yml --core 50
$snakemake -s classifier_quant.smk --configfile classifier_quant.yml --core 50
```



Author: LINA LU	20220815

