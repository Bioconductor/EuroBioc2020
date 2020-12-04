---
layout: default
---

{% include header.md %}

# Birds Of a Feather

## Infrastructure to enable spatially resolved transcriptomics data analysis with Bioconductor

_Organizers_: Dario Righelli, Helena L. Crowell, Lukas M Weber

_Abstract_: 
Spatially resolved transcriptomics technologies allow joint collection of transcriptome-wide mRNA expression and spatial coordinates on a tissue slide. Multiple technologies are emerging (including 10x Genomics Visium, seqFISH, and MERFISH), which has led to several Bioconductor contributors actively developing object classes and methods for data analysis. In light of this, we propose an open discussion on what has been developed so far to provide a state-of-the-art analysis workflow, while delineating future directions.

Goal 1: SpatialExperiment object class and infrastructure

The primary goal for this meeting will be to consolidate ideas regarding how to implement flexible data object classes and infrastructure to accommodate various types of spatially resolved transcriptomics data. This will be based around the existing SpatialExperiment Bioconductor classes (https://bioconductor.org/packages/SpatialExperiment). In particular, there is a need to support additional features in the SpatialExperiment classes, including: 

Multiple samples accommodation.
Storage of images for multiple samples and of different resolutions.
Allow images to be stored as paths and/or URLs (for later retrieval) and/or shareable data objects (e.g. grobs) for direct accession.

Outcome: After agreeing on a class design, modular tasks could be collected as GitHub issues, taken up by different developers at a later point in time, and merged through pull requests.
In addition, we aim to decide on a strategy to draft a paper describing SpatialExperiment. For example, this could take the form of a short paper describing the class structure and showcasing examples from several technologies. This could be drafted as a collaborative paper, with sections assigned to developers working with each technology (e.g. Visium, seqFISH).

Goal 2: Case-study analysis workflows 

A secondary goal will be to develop vignette-style analysis workflows to be built around the SpatialExperiment object classes. These could be integrated within a larger online textbook on spatial resolved transcriptomics analysis with Bioconductor (currently under development at https://github.com/lmweber/OSTA-base). Examples include:

Single-sample scenario where it is of interest to carry out an in-depth analysis
including, e.g., basic preprocessing, clustering and cell-subpopulation assignment and/or deconvolution using a reference dataset with single-cell resolution, identifying spatially variable features across the tissue and/or within each subpopulation, etc.
Multi-sample scenario, e.g., where biological replicates from different experimental conditions or samples from different tissues are available.
Examples using different technologies (e.g. Visium, seqFISH, etc).

Outcome: Identify missing workflow examples in the current materials, and assign to developers working with each technology, to add as new workflow chapters in OSTA or other resources (e.g. GitHub repositories). 

While the conversation will be kept in an open format to enable participation from attendees coming from diverse backgrounds and experience levels, as already mentioned, we would like to document and structure the output of this BoF as one or more documents to provide headlines and points of view on spatially resolved transcriptomics data analysis in Bioconductor.

## Creating Github Actions for building and checking Bioconductor Packages

_Organizers_: Mike L Smith, Constantin Ahlmann-Eltze

_Abstract_:
Many developers use GitHub Actions to build and test their packages outside of the Bioconductor build system.  Although there is an existing ecosystem of R related GitHub actions, as far as we are aware there are no Bioconductor specific actions or documentation.

In this birds-of-a-feather session we hope to discuss whether developing these is something that would benefit the community, and what would be required to make a useful resource for package developers.

## The emerging R ecosystem for microbiome research

_Organizers_: Felix G.M. Ernst, Sudarshan A Shetty, Ruizhu Huang, Domenick James Braccia, Hector Bravo, Leo M Lahti

_Abstract_:
Dedicated class structures are fundamental for the R/Bioconductor ecosystem. In microbiome research, the phyloseq class has become a commonly accepted standard. However, recent developments in Bioconductor classes and growing sample sizes have opened up new opportunities and requirements to extend this popular format in order to address emerging research needs. Extended support is needed for instance for hierarchically structured data, spatio-temporal variation, linking with other data types, and for general performance optimization.

This birds-of-a-feather session aims to map and discuss the current state of the microbiome research ecosystem in R/Bioconductor. In addition, we will discuss the possibilities to extend the current state-of-the-art by taking advantage of the TreeSummarizedExperiment and SingleCellExperiment classes, which provide thoroughly tested tools for hierarchical data, spanning both sample and feature spaces. This can bring extended support for topics such as sparse matrices and multiple assays while providing improvements in speed and memory compared to the currently available solutions for microbiome data. These extensions would be compatible with the widely used phyloseq class as well as other raw data types, thus enabling seamless conversion. By linking microbiome data more tightly to other already established Bioconductor classes, we hope to reduce overlapping development efforts, improve the interoperability of available tools, and ensure the long-term sustainability of the ecosystem. The class structure and associated packages are now under active development as a joint effort of multiple research teams. Feedback and contributions from the broader community will be highly valued.

The event consists of an introduction followed by facilitated discussion. The event may be splitted in smaller sub-teams  depending on the number of participants. The planned outputs of the session include 1) summarizing the requirements for a sustainable microbiome data structure; 2) mapping of the current and emerging needs for data analytical tools; and 3) a roadmap for collaborative development.

## Ten simple rules for thriving in bioinformatics research

_Organizers_: Federico Marini, Charlotte Soneson, Davide Risso

How is it possible as bioinformaticians/computational biologists to develop our own career in the Life Sciences, e.g., by developing software?

Focusing on the fields of academia, consulting, and research service, we initiated a discussion at BioC2020, where many participants have provided valuable inputs that nicely steered the conversation.

Now, we would like to summarize these experiences and suggestions and distill them into a final tentative list of “Ten simple rules” (a well-known and insightful format for the readers of Plos Computational Biology), potentially responding to https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007531, where we’d like to highlight the research aspects over the pure “research support”.

We’d like to discuss opportunities for enhancing recognition of our work (e.g. by being able to publish workflows and software); encourage evaluation based on additional metrics than the “traditional” papers and citations; how can we acquire funding being mostly in the software/method development. We welcome contributions from attendees of all possible backgrounds and academic levels (industrial/academic, students/PIs).

The desired outcome would be structured in a Google Slideset, for collaborative live editing, where we aim to define the 10 directions to touch upon when writing the accompanying manuscript (for which a skeleton is set up with the Manubot framework, https://github.com/drisso/ThrivingInBioinformatics).

# Workshops

## tidytranscriptomics : introduction to tidy analysis of single-cell and bulk RNA sequencing data

_Duration_: 90 mins

_Authors_: Stefano Mangiola, Maria Doyle

_Abstract_: 
This workshop will present how to perform analysis of RNA sequencing data following the tidy data paradigm. The tidy data paradigm provides a standard way to organise data values within a dataset, where each variable is a column, each observation is a row, and data is manipulated using an easy-to-understand vocabulary. Most importantly, the data structure remains consistent across manipulation and analysis functions.

This can be achieved for RNA sequencing data with the tidySCE, tidybulk, tidyHeatmap and tidyverse packages. The tidybulk package provides a tidy data structure and a modular framework for bulk transcriptional analyses, the tidySCE package provides similar for single-cell transcriptional analyses, and tidyHeatmap provides a tidy implementation of ComplexHeatmap. These packages are part of the tidytranscriptomics suite that introduces a tidy approach to RNA sequencing data.

## Quantitiative analysis of ChiP-seq, ATAC-seq, and related DNA enrichment assays

_Duration_: 90 mins

_Author_: Rory Stark

_Abstract_:
This workshop will demonstrate the steps involved in performing a 
quantitative analysis of ChIP-seq data in _Bioconductor_
(up to and including differential binding analysis), 
with some discussion of related assays such as ATAC-seq.
Particular attention will be paid to processing of aligned reads, 
including blacklisting, greylisting, filtering for quality and duplication,
and the particular challenges presented when normalizing these data.
While the workshop follows the _DiffBind_ package vignette, the use of 
a number of other _Bioconductor_ packages is discussed, including
_csaw_, _ChIPQC_, _edgeR_, _DESeq2_, and 
_GreyListChIP_.

## Integration of ChIP-seq and RNA-seq data in R

_Duration_: 45 mins

_Authors_: Mahmoud Ahmed

_Abstract_:
Researchers use ChIP binding data to identify potential transcription factor binding sites. Similarly, they use gene expression data from sequencing or microarrays to quantify the effect of the factor overexpression or knockdown on its targets. The integration of the binding and expression data therefore can be used to improve the understanding of a transcription factor function. In this workshop, we present a complete workflow for integrating the gene expression (RNA-seq) and DNA-binding data (ChIP-seq) to predict the combined function of two transcription factors using R/Bioconductor. The example we will be using in the workshop is from real datasets of two functionally and evolutionary related transcription factors YY1 and YY2 in HeLa cells. We will try to identify the factor-specific and the shared targets of the factors in this particular cell line. Then we will use a technique to find out the aggregate functions of the factors on their individual (inducer or repressor) and common targets (cooperative or competitive). The first half of the workshop would be dedicated to introduce and walk through the workflow interactively in a live demo. In the second half, we will explore the output of the analysis to answer specific questions in a quiz format followed by an open discussion. 

## NewWave, new R package for dimensional reduction and batch effect removal for single cell RNA-seq data

_Duration_: 45 mins

_Authors_: Federico Agostinis, Chiara Romualdi, Gabriele Sales, Davide Risso

_Abstract_:
The fast development of single cell sequencing technologies in the recent years has generated a gap between the throughput of the experiments and the capability of analizing the generated data. 
In this package, NewWave,  we implement mini-batch stochastic gradient descent and the possibility to work with HDF5 files. We decided to use a negative binomial model following the observation that droplet sequencing technologies do not induce zero inflation in the data. Thanks to these improvements and the possibility of massively parallelize the estimation process using PSOCK clusters, we are able to speed up the computations with the same or even better results than zinbwave. This type of parallelization can be used on multiple hardware setups, ranging from simple laptops to dedicated server clusters. This, paired with the ability to work with out-of-memory data, enables us to analyze datasets with milions of cells.

## preciseTAD: a machine-learning framework for predicting boundaries of 3D genomic elements

_Duration_: 45 mins

_Authors_: Spiro C Stilianoudakis, Mikhail Dozmorov

_Abstract_:
High-throughput chromosome conformation capture technology (Hi-C) revealed extensive DNA folding into discrete 3D structures referred to as Topologically Associating Domains (TADs) and loops.  TADs are critical for cellular processes like gene regulation and cell differentiation. The relatively low resolution of Hi-C data (tens of kilobases in size) prevents precise mapping of TAD boundaries by conventional TAD-callers. In contrast, the enrichment of high-resolution genomic annotations at boundaries (ChIP-seq, ~100 bases) offers a computational approach toward improved boundary identification. The workshop describes how to transform TAD-calling into a supervised machine learning framework. We demonstrate how to build and optimize a random forest model that prioritizes known molecular signatures of domain boundaries to predict each base's probability of being a boundary. We utilize density-based clustering (DBSCAN) and partitioning around medoids (PAM) to identify boundary regions and the most probable points. Boundaries identified by preciseTAD show strong signal enrichment of known boundary signatures compared with Arrowhead- and Peakachu boundaries. This workshop uses the R programming environment, the preciseTAD R package (https://bioconductor.org/packages/preciseTAD/), and can be performed with most operating systems on a single computer.

## Spatial Transcriptomics workflow using SpatialExperiment package classes

_Duration_: 45 mins

_Authors_: Dario Righelli, Lukas M Weber, Helena Lucia Crowell, Davide Risso

_Abstract_:
Cells spatial organization within tissues is important for the understanding of cell communications and interactions in developmental, physiological, and pathological states.
Lately several technologies are emerging for the joint extraction of transcripts quantification and their spatial organization at single cell level.
Among others, seqFISH and 10x Visium Spatial Gene Expression (VSGE) seem to be the most widely used platforms because of their ability to extract hundreds to thousands of transcripts at the same time in the first case and for the standard protocol implemented in the latter.
Briefly, seqFISH combines temporal barcodes in multiple hybridization rounds with microscopy fluorescent imaging to detect the spatial coordinates of transcripts.
On the other hand, VSGE allows the joint analysis of microscopy images of a chip-placed frozen piece of tissue and its transcriptome at spot-barcode level.
Given the lack of dedicated general tools for spatial transcriptomics, we developed the SpatialExperiment Bioconductor package, which presents a set of two dedicated classes (SpatialExperiment and VisiumExperiment) both extending the SingleCellExperiment Bioconductor class with slots and methods for spatial transcriptomics data handling.
This work could become the basis for facilitating further development of methods aimed to analyze this type of data, e.g., for spatial clustering and the identification of spatially variable genes.
Here, we present a short workshop showing how to analyze spatial transcriptomics with our SpatialExperiment package and how to use it in combination with other analysis software.
We will start by showing how to combine a SpatialExperiment with a SingleCellExperiment to store multimodal datasets such as seqFISH and scRNAseq and how to manage them to graphically visualize the data.
Additionally, we will show how to create a VisiumExperiment object starting from 10x Genomics public datasets and how to manage it for clustering annotation and plotting.

## SPEAQeasy: a Scalable Pipeline for Expression Analysis and Quantification for R/Bioconductor-powered RNA-seq analyses

_Duration_: 90 mins

_Authors_: Nicholas J Eagles, Emily E Burke, Jacob Leonard, Brianna K. Barry, Joshua M. Stolz, Louise Huuki, BaDoi N. Phan, Violeta Larios Serrato, Everardo Gutiérrez-Millán, Israel Aguilar-Ordoñez, Andrew E. Jaffe, Leonardo Collado Torres

_Abstract_:
RNA sequencing (RNA-seq) involves a large number of individual steps before yielding directly valuable information, such as differential gene expression data, from raw reads. Existing software tools are typically specialized, only performing one step-- such as alignment of reads to a reference genome-- of a larger workflow. The demand for a more comprehensive and reproducible workflow has led to the production of a number of publicly available RNA-seq pipelines. However, we have found that most require computational expertise to set up or share among several users, are not actively maintained, or lack features we have found to be important in our own analyses. In response to these concerns, we have developed a Scalable Pipeline for Expression Analysis and Quantification (SPEAQeasy), which is easy to install and share, and provides a bridge towards R/Bioconductor downstream analysis solutions. SPEAQeasy is user-friendly and lowers the computational-domain entry barrier for biologists and clinicians to RNA-seq data processing as the main input file is a table with sample names and their corresponding FASTQ files. SPEAQeasy is portable across computational frameworks (SGE, SLURM, local, docker integration) and different configuration files are provided.

In this workshop, we would first provide an overview of what steps SPEAQeasy performs. Participants would then be guided through a complete end-to-end RNAseq analysis, leveraging SPEAQeasy. This would involve downloading publicly available example data, configuring SPEAQeasy for the experiment, and running it. Afterward, we would demonstrate how SPEAQeasy outputs can be coupled with user-provided genotype calls to resolve identity issues, perform a differential expression analysis, and examine gene ontology results using Bioconductor packages (SummarizedExperiment, limma, clusterProfiler, among others). The goal would be to introduce a straightforward but flexible workflow for researchers at any level of technical experience.

## iSEE: Interactive SummarizedExperiment Explorer

_Duration_: 90 mins

_Authors_: Kevin Rue, Charlotte Soneson, Federico Marini

_Abstract_:
The iSEE (Interactive SummarizedExperiment Explorer) Bioconductor software package provides a general visual interface for exploring data in a SummarizedExperiment object, using the shiny R package.

Following a brief introduction on SingleCellExperiment objects and an overview of the iSEE software functionality, this workshop is structured in self-contained recipes and demonstrates via simple hypothesis-driven tasks how to use iSEE to create, configure, and interact with the user interface, for the interactive exploration and interpretation of biological data sets to answer scientific hypotheses and generate insights into biological phenomena.

Each recipe shows the final result of the target iSEE configuration, and provides simple hints to reach the desired outcome, as well as full descriptions of the solutions, both by pointing-and-clicking in the web application (for interactive configuration on the fly), and by programmatic approaches on the console (using simple code scripts for preconfigured setup).

As a complement to the hands-on session, instructors will be available to answer questions facilitating the usage and configuration of the iSEE package for the predefined exercises, with some scope for participants to apply iSEE to their own data sets, and fuel the discussion with more questions about specific use cases.

## QFeatures: Quantitative features for mass spectrometry data

_Duration_: 45 mins

_Authors_: Laurent Gatto, Christophe Vanderaa

_Abstract_:
The QFeatures package provides infrastructure (that is classes to store data and the methods to process and manipulate them) to manage and analyse quantitative features from mass spectrometry experiments. It is based on the SummarizedExperiment and MultiAssayExperiment Bioconductor classes. Assays in a QFeatures object have a hierarchical relation: proteins are composed of peptides, themselves produced by spectra. Throughout the aggregation and processing of these data, the relations between assays are tracked and recorded, thus allowing users to easily navigate across spectra, peptide and protein quantitative data. In this workshop, we will demonstrate how to import data as QFeatures objects, how to process and analyse data in QFeatures objects, and how to interpret the results. Some familiarity with Bioconductor data analysis, in particular the SummarizedExperiment class, is recommended to follow this short workshop. 


## Multi-omic Integration and Analysis of cBioPortal and TCGA data with MultiAssayExperiment

_Duration_: 45 mins

_Authors_: Marcel Ramos, Levi Waldron

_Abstract_:
This workshop demonstrates the leveraging of public multi-omics databases, such as cBioPortal and The Cancer Genome Atlas (TCGA), through the use of the cBioPortalData and curatedTCGAData experiment data packages. It provides users with the basics of data management, using the MultiAssayExperiment data class and the TCGAutils utility package, and example analyses of multiple assays associated with a single set of biological specimens. In addition to providing a basic overview of key data classes, such as MultiAssayExperiment and RaggedExperiment, this workshop intends to provide an overview of cBioPortalData and curatedTCGAData experiment data packages and TCGAutils functionality aimed at enhancing the ease-of-use of TCGA data.

## Seamless Integration of Mass Spectrometry Data from Different Sources with the Spectra Package

_Duration_: 45 mins

_Authors_: Johannes Rainer, Sebastian Gibb, Laurent Gatto

_Abstract_:
Mass spectrometry (MS) data is a key technology in modern proteomics and metabolomics experiments. Due to continuous improvements in MS instrumentation, the generated data can easily become very large. Also, different additional resources of MS data exist, such as spectra libraries and databases, all with their own specific file formats that sometimes do not support manipulations of the original data.

Learning from experiences with the MSnbase Bioconductor package we developed a novel infrastructure to handle MS spectral data in R, the Spectra package. This package implements a clear separation of user functionality from code to provide, store and import MS data. Different backends can hence be used that enable access to data from various data resources or that are designed specifically for very large MS data sets. Data manipulations are by default not directly applied to the data but cached in a lazy processing queue which allows analyses also of read-only data representations.

This workshop shows the expandability of the new infrastructure to enable a seamless integration and analysis of MS data from a variety of input formats illustrated by a simple matching of experimental MS2 spectra against a public spectral database and export of the data in a format commonly used for exchange of MS2 data.

## Trajectory inference across conditions: differential expression and differential progression

_Duration_: 90 mins

_Authors_: Hector Roux de Bezieux, Kelly Street, Koen Van den Berge

_Abstract_:
In single-cell RNA-sequencing (scRNA-seq), gene expression is assessed at the level of single cells. In dynamic biological systems, it may not be appropriate to assign cells to discrete groups, but rather a continuum of cell states may be observed, e.g. the differentiation of a stem cell population into mature cell types. This is often represented as a trajectory in a reduced dimension of the scRNA-seq dataset.

Many methods have been suggested for trajectory inference. However, in this setting, it is often unclear how one should handle multiple biological groups or conditions, e.g. constructing and comparing the differentiation trajectory of a wild type versus a knock-out stem cell population.

In this workshop, we will explore methods for comparing multiple conditions in a trajectory inference analysis. We start by integrating datasets from multiple conditions into a single trajectory. By comparing the conditions along the trajectory’s path, we can detect large-scale changes, indicative of differential progression. We also demonstrate how to detect subtler changes by finding genes that exhibit different behaviors between these conditions along a differentiation path.

# Long talks (20 minutes)

## distinct: a novel approach to differential analyses

_Speaker_: Simone Tiberi

_Abstract_:
We present distinct, a statistical method to perform, via hierarchical permutation tests, differential analyses between groups of densities. distinct is a general and flexible tool: due to its fully non-parametric nature, which makes no assumptions on how the data was generated, it can be applied to a variety of datasets. It is particularly suitable to perform differential analyses on single cell data, such as single cell RNA sequencing (scRNA-seq) and high-dimensional flow or mass cytometry (HDCyto) data. While most methods for differential expression target differences in the mean abundance between conditions, single-cell data can show more complex variations. distinct, by comparing full distributions, identifies, both, differential patterns involving changes in the mean, as well as more subtle variations that do not involve the mean (e.g., unimodal vs. bi-modal distributions with the same mean). distinct explicitly models the variability between samples (i.e., biological replicates), can adjust for covariates (e.g., batch effects) and allows multi-group (i.e., >2 groups) comparisons. We will present results, based on scRNA-seq and HDCyto simulated and experimental datasets, where distinct outperforms several competitors and is able to detect more patterns of differential expression compared to canonical differential methods. distinct is freely available as a Bioconductor R package.

## pARI package: valid double-dipping via permutation-based All Resolutions Inference

_Speaker_: Angela Andreella

_Abstract_:
The cluster extent-based thresholding (Woo et al., 2014) for functional Magnetic Resonance Imaging (fMRI) data is widely popular for finding neural activation associated with some stimulus. Contiguous collection of voxels, i.e., cluster, is declared significant at a pre-specified level; however, it suffers from the spatial specificity paradox (Woo et al., 2014). Since the method tests the hypothesis that none of the voxels in the cluster are active, rejecting this null hypothesis only allows the claim that there is at least one active voxel inside the cluster. The number of active voxels and their spatial location remains unknown, and making follow-up inference inside the cluster (``drilling down) leads to a double-dipping problem and inflated Type I error rate. For that, Rosenblatt et al. (2018) develop a closed testing method (Marcus et al., 1976) with Simes local tests (Sarkar, 2008) called All-Resolution Inference (ARI), to infer the number of true discoveries simultaneously for all possible subsets of hypothesis still controlling the familywise error rate (FWER). The method associates for each cluster the lower bound of the number of truly active voxels. However, ARI can lose power if the data are strongly correlated, e.g., fMRI data or genetic data. For that, we re-phrase this method using the permutation theory. Permutation theory adapts to the spatial correlation structure, gaining power over parametric approaches (Pesarin and Salmaso, 2010). Therefore, we present pARI: an R package developed to compute the permutation-based All-Resolution Inference (ARI). The main function, called pARIbrain, takes as input a list of copes, i.e., contrast maps, one for each subject, given by neuroimaging tools as FSL, SPM, etcetera. The user can then insert a cluster map from some neuroimaging software, a region of interest (ROI) image, or any subsets of voxels of interest. It is possible to construct these cluster maps using a supra-threshold statistic rule; the user can specify the threshold into the arguments of the function pARIbrain. Therefore, the function pARIbrain returns the lower bounds of true discoveries, i.e., active voxels, for each cluster specified. The package was developed mainly for the fMRI scenario; however, we develop also the function pARI. It takes as inputs the permutation null distribution and the indexes of the hypothesis of interest . pARI then returns the lower bound for the number of true discoveries inside the set of hypotheses specified. Besides, thanks to the permTest and signTest function, the user can easily compute the permutation null distribution concerning the two-sample t-tests and one-sample t-tests. The main functions are written in C++ to improve the computational efforts associated with the computation of the permutation p-values null distribution. We applied the method in fMRI, electroencephalogram (EEG), and genetic data. The package is available on https://github.com/angeella/pARI.

## Large-Scale Topological Changes Restrain Malignant Progression in Colorectal Cancer

_Speaker_: Alejandro Reyes

_Abstract_:
When diagnosing cancer, pathologists use nuclear morphology as a hallmark of tumor cells, but the topological changes linked to alterations in the shape of the nucleus were poorly understood. To uncover the molecular changes that occur in the 3D structure of cancer genomes, we integrated topological maps for colon tumors and normal colon tissues with epigenetic, transcriptional, and imaging data. We characterized tumor-associated changes in chromatin loops, topologically associated domains (TADs), and large-scale compartments. We found that chromatin loop rewiring contributes to oncogenic gene expression programs. Whereas TAD structures are largely stable, the spatial partitioning of the open and closed genome compartments is compromised in tumors and is accompanied by hypomethylation and histone mark rearrangements. We also identified a compartment that is reorganized in tumors and is at the interface between the canonical A and B compartments. Similar changes were evident in non-malignant cells that had accumulated excessive cellular divisions. Our analyses suggest that these compartment changes repress stemness and invasion programs while inducing immunity genes, and they may therefore restrain malignant progression. By shedding light on a protective mechanism that acts through 3D genome organization, the findings can potentially be used to explore new cancer therapies.

## scGCN: A Geometric Deep Learning Framework on Single-cell Gene Networks

_Speaker_: Elyas Heidari

_Abstract_:
Motivation: Single cell transcriptomics have enabled resolution of transcriptional diversities of individual cells, hence the possibility to classify cells’ identity and function via their transcriptional state. Current single-cell RNA sequencing (scRNA-seq) classification methods use genes expression as independent features to assign cell types, disregarding their interactome which is of great importance in consistency of different cell states. Cell type classification based on a few marker genes, for example, yields better results in intra data set settings where both train and test cells are sampled from the same dataset (Abdelaal et al. 2019) and may fail to correctly classify cell types under perturbations when expression of the predefined marker genes varies with respect to the control reference data set. Moreover, using genes expression as independent features can delude cell lineage relationships for cell states and cell types with subtle differences in their gene expression profiles. Method: We propose, scGCN, a deep learning framework for analysis of single-cell RNA-seq data. Our framework consists of three steps; 1) adaptation of a geometric deep learning model applicable to gene regulatory networks (GRNs), 2) construction of a modular GRN, which serves as the input graph structure for the geometric deep learning model, and 3) validation of the model on available scRNA-seq data sets. For each of these steps, we face a multitude of possible methods. We would like to systematically explore the solution space in order to find the best amongst the set of possible solutions. This brought us to first integrating all steps into one pipeline, and then optimizing the whole pipeline to achieve the best solution. For instance, for GRN reconstruction, we employ probabilistic graphical models, a priori known transcription factor interactions, or a combination of both. Results: We employed two pairs of scRNA-seq datasets from two organisms, that is, a pair of distinct datasets on human peripheral blood mononuclear cells (PBMC) and another pair on Zebrafish embryo, we demonstrate the potentials of scGCN for transfer learning. We came up with a parameter path consisting of 1) selection of highly variable genes, 2) reconstruction of GRN with Gaussian Graphical Models (partial correlation networks), 3) GRN classification with InceptionGCN, suitable for transferable learning on scRNA-seq data, as we tried it on two datasets. Our results show as high but more stable accuracy match between the test and train datasets compared to a fully connected network. Conclusion: We illustrate a two-fold benefit of graph representation of cell states in removing unwanted and technical variations (i.e., batch effects) as well as the potentials for transferable learning, i.e., learning models which are least affected by data set specific variations and hence are more general and can be used for predictions on new data. Whereas unwanted variations among different data sets makes any supervised learning method prone to overfitting and dataset specific model artifacts, we suggest that a supervised learning approach on GRNs imposes the proper inductive bias (Baxter 2000, Neyshabur, Tomioka, and Srebro 2014) for ‘transferable learning’ (Levie, Isufi, and Kutyniok 2019), such that an available data set can be used to train a model which is able to predict cell types in new (unseen) data sets later on. To the best of our knowledge, our proposed method is the first end-to-end pipeline for network construction and graph deep learning, all together. Moreover, this is the first pipeline employing graph neural networks in the field of single-cell transcriptomics to perform cell-type annotation. Availability: We implement our pipeline as a software package in R and python. The package `scGCNUtils` provides the user with the building blocks of the pipeline as well as performing pre-/post- exploratory data analysis (EDA). The user is able to run an end-to-end pipeline by defining their desired parameters, and decide on the ultimate parameter path using EDA tools provided for each module. The package also facilitates benchmarking, in that each step of the pipeline is implemented as an independent module and can be tweaked or substituted by another module with the same functionality. In the end, the performance and efficiency of the path can be validated by the package. Logging and caching are added to the package for the same purpose. The package is designed to meet state-of-the-art reproducibility standards and criteria of open science, such as coherence, integrity, documentation, readability, and testability. The package is going to be submitted to Bioconductor shortly. The development version of code is available at: https://github.com/EliHei2/scGCN

## Infinity Flow: High-throughput single-cell quantification of 100s of proteins using conventional flow cytometry and machine learning

_Speaker_: Etienne Becht

_Abstract_:
Objective: Systems immunology requires tools that can reveal the complex interplay of a wide array of cell types in health and disease. Modern flow cytometers can quantify around 20 proteins at the single-cell level, allowing a large number of cell clusters to be identified in tissues. However, it remains challenging to exhaustively map each of these clusters to known cell types. In addition, for antibody panel designing, inclusion of well-characterized markers conflicts with inclusion of more exploratory markers. Higher-dimensional alternatives to flow cytometry exist, such as mass cytometry, or oligonucleotides-tagged antibodies - but they are limited in terms of cellular throughput, higher in cost and not as widely available. To overcome these limitations, we propose a new computational workflow to fully leverage massively-parallel cytometry experiments. Methods: Experimentally, cells are stained with a customizable 10-20 antibodies “backbone” panel which defines the cellular population structure of the sample. This is followed by aliquoting into hundreds (200-300) of individual wells, each containing a unique “exploratory” antibody. Data acquisition is performed using conventional flow cytometry. Computationally, multivariate non-linear regression models (neural networks, boosted trees or support vector machines) are trained, one for each well. These models learn how to impute an exploratory antibody signal from backbone measurements. The models are each applied to the backbone data to obtain hundreds of imputed measurements across millions of single cells. Results: We applied this workflow to single-cells dissociated from the lungs of C57/B6 mice, resulting in a data matrix of 2,660,000 events with imputed expression levels for 266 proteins. Imputed expression and co-expression patterns are accurate and consensual across regression models, allowing near-exhaustive annotation of both immune and non-immune cell types (32 out of 33 clusters annotated, 96.7% of cells), and finer granularity as compared to dimensionality reduction or clustering. Further, dimensionality reduction reveals additional biologically meaningful heterogeneity when performed on imputed data compared to the backbone data. Conclusion: By enabling broad characterization of cellular networks at the protein level, this pipeline is well-suited for systems immunology approaches. Given the wealth of data generated by this pipeline and its accessibility, we anticipate it will be widely adopted. The computational pipeline was recently published as a Bioconductor package allowing easy reproduction of this approach on any suitable dataset.

## Evaluation of cell populations in the locus coeruleus through spatial transcriptomics

_Speaker_: Abby Spangler

_Abstract_:
The locus coeruleus (LC) is the main site for the synthesis of norepinephrine in the human brain and is therefore involved in the physiological responses to stress and panic. It is indicated in neurophysiological disorders such as Parkinson’s, Alzheimer’s, clinical depression, panic disorder and anxiety. It is composed mostly of medium-sized neurons however the complete spatial diversity of cell-types and gene expression that compose this area of the brain is still unknown. In this study, we perform spatial RNA-seq using the 10x Genomics Visium platform on two slices of LC dissected from fresh frozen human brain. The two slices are from the same donor and are considered spatial replicates taken 100uM apart. We use the spatialLIBD and other Bioconductor packages to cluster and visualize the data in order to assess spatial patterns of gene expression and new cell populations, which are relatively unknown compared to the dorso lateral prefrontal cortex (Maynard et al, bioRxiv, 2020), thus making the LC a more challenging section of the brain to study.

## Standardised and reproducible analysis of mass spectrometry-based single-cell proteomics data

_Speaker_: Christophe Vanderaa

_Abstract_:
Recent advances in sample preparation, processing and mass spectrometry (MS) have enabled the emergence of MS-based single-cell proteomics (SCP). However, the lack of computational framework limits the optimization and wider application of these new technologies. We aim to fill this gap and propose a standardized pipeline for the analysis of MS-based SCP data. Our packages rely on existing Bioconductor infrastructure, such as classes and functions for single-cell RNA sequencing and MS-based quantitative proteomics. The first package, scpdata, disseminates curated SCP data sets for method development and benchmarking. The second package, scp, implements functions to streamline the analysis of SCP data. We demonstrate the application of our pipeline by replicating and improving on the analyses of two published data sets and show its relevance in the processing and interpretation of MS-based SCP data.

## The effect of phosphorylation on the spatial proteome

_Speaker_: Oliver Crook

_Abstract_:
 For proteins to function correctly they must be localised to the correct subcellular compartment. Phosphorylation is a molecular toggle that can alter the function of a protein. However, this phenomenon is not well understood. We present a high-throughtput mass-spectrometry method to obtain simultaneous localisation information of phosphorylated and non-phosphorylated proteins. Using a Bayesian toolkit, we determine whether modification has an effect on protein subcellular localisation. Our results provide insights into the functional role of phosphorylation.

# Short talks (8 minutes)

## Doublet identification and characterization in single-cell data with scDblFinder

_Speaker_: Pierre-Luc Germain

_Abstract_:
Doublets (i.e. two cells captured in the same reaction volume and hence appearing as a single cell) can form at a high frequency in single-cell sequencing studies (often 10-20%). They can appear as spurious 'new cell types' or distort trajectory and co-expression analyses. While experimental strategies can be used to mitigate this effect, they are insufficient and need to be complemented with doublet identification methods. In this talk I discuss the main doublet identification strategies, and present the scDblFinder package, its main approach, and some lessons learnt from its implementation. I then demonstrate its superiority (in terms of both speed and accuracy) to alternatives on a vast majority of datasets, and finally turn to the problem of doublet characterization and doublet enrichment analysis.

## A new statistical approach for the simultaneous clustering of genes and cells in spatial transcriptomic experiments

_Speaker_: Andrea Sottosanti

_Abstract_:
In the last years, we witnessed a substantial improvement in the efficiency of DNA sequencing technologies with the rise of new advanced protocols for spatial transcriptomics, such as seqFISH and 10X-Visium. With respect to the classical transcriptomic technologies, these tools can provide, starting from a tissue sample, a full reconstruction of the spatial location of the cells, and the expression of thousands of genes inside each cell. The increasing popularity of such advanced technologies has motivated the development of statistical methods that aim to identify and study the so-called spatially expressed genes, i.e. genes whose expression in a cell affects the expression in the surrounding ones. These kinds of genes have a great scientific interest, as comprehending their functions and their interactions across multiple cells might lead to a deeper understanding of several complex biological mechanisms. Another relevant aspect of scientific interest in the analysis of a tissue is the classification of the cells. Distinguishing, for example, a tumor cell from a stromal or an immune cell is vital. Again, statistics, and in particular clustering techniques, plays a core role in this field. Including the additional spatial information of the cells carried by the new transcriptomic technologies would mean a big step forward in the concrete ability of a classifier of distinguishing groups of cells of different nature. In this talk, we present a new advanced statistical method for the analysis of complex genomic datasets provided by the 10X-Visium technology. We investigate the properties of the matrix variate distributions, a simple though theoretically coherent way to define a statistical model over the entire data matrix, capturing the complex dependency structure of the data. We further assume that both genes and cells can be clustered into groups, dividing the data matrix into rectangular, non-overlapped blocks. This procedure, known is statistical literature as co-clustering, allows us to classify the nature of cells and to detect groups of spatially expressed genes only in some specific cell types.

## Integration of preprocessed single-cell datasets with bulk differential analysis results for DESeq2

_Speaker_: Kwame Forbes

_Abstract_:
The Bioconductor project contains numerous single-cell RNA-seq datasets readily available for download, which can provide insight into the expression of genes in particular groups of cells, for example cells grouped by cell type. Here, we have developed a simple interface for users with bulk RNA-seq differential expression results tables to make use of public single-cell datasets, by visualizing the cell-type specific expression patterns of their genes of interest. This new interface leverages metadata-gathering functionality in the tximeta Bioconductor package to automatically map across transcript or gene identifiers. We also implement functionality to determine if significant bulk gene sets tend to be expressed in specific groups of cells in the single-cell data, via visualization and contingency table-based hypothesis testing. The new functions are incorporated into the devel branch of the DESeq2 package and will be released in October 2020.

## Using VisiumExperiment objects at spatialLIBD package

_Speaker_: Brenda Pardo

_Abstract_:
spatialLIBD is a package to interactively visualize the LIBD human dorsolateral pre-frontal cortex (DLPFC) spatial transcriptomics data (Maynard, Collado-Torres, Weber, Uytingco, et al., 2020). It contains functions for (1) accessing the spatial transcriptomics data from the LIBD Human Pilot project generated with the Visium platform from 10x Genomics, (2) visualizing the spot-level spatial gene expression data and clusters, and (3) inspecting the data interactively either on the user’s computer or through a web application. spatialLIBD used to employ R objects from the SingleCellExperiment class to store the data, nevertheless, Righelli et al created a more accurate class, VisiumExperiment, for spatial transcriptomics data. In this work, we made spatialLIBD able to use objects from the VisiumExperiment class by constructing a multi-sample object and modifying the spatialLIBD functions to support VisiumExperiment objects as an input.

## Normalizing ChIP and ATAC sequencing data for differential analysis

_Speaker_: Rory Stark

_Abstract_:
Bioconductor includes a number of packages for performing differential analysis of ChIP-seq and ATAC-seq DNA enrichment assays. The most referenced differential binding analysis packages, DiffBind and csaw, rely on underlying differential modelling packages designed for RNA-seq (edgeR and DESeq2). However, assumptions made when normalizing RNA-seq data, notably that there are a core set of sites that do not change binding affinity, are generally not applicable to DNA enrichment data. As a result, the normalization step plays a crucial role in identifying differentially bound sites that reflect the underlying biology. In a peak-based analysis, using the native normalization methods in the RNA-seq modelling packages directly on a consensus read matrix will result in conclusions contrary to the biology being studied if there are systematic shifts in transcription factor binding or chromatin structure. We undertook a study of the impact of the more than a dozen different normalization schemes now possible within DiffBind, including background binning, spike-ins, and parallel factors. This study demonstrates that the reference set of reads used as the basis for normalization is more important than the specific normalization algorithm (e.g. TMM or RLE) or differential modelling package (edgeR or DESeq2). In addition to normalization, the use of blacklists, greylists, and the handling of duplicate reads can also have a major impact on differential results; given a “long” talk, I will discuss each of these steps, using a novel dataset that includes UMIs to identify “true” duplicates.

## Methrix: a package for systematic aggregation and analysis of bisulfite sequencing data

_Speaker_: Anand Mayakonda

_Abstract_:
Whole genome bisulfite sequencing (WGBS), measures DNA methylation at base-pair resolution resulting in large bedGraph like coverage files. Fundamental downstream analyses require summarization of coverage files into methylation and coverage matrices, whose dimensions rapidly increase along with the number of sequenced samples. Current options for processing such files are hindered by discrepancies in the file format specification, speed, and memory requirements. We developed methrix, an R package, which provides a toolset for systematic analysis of large datasets. Core functionality of the package includes a comprehensive bedGraph or similar tab-separated text file reader - which summarizes methylation calls based on annotated reference indices, infers and collapses strands, and handles uncovered reference CpG sites while facilitating a flexible input file format specification. Additional optimized functions for quality control filtering, sub-setting, and visualization allow user-friendly and effective processing of WGBS results. Additional arguments are available for working with large (n>100) cohorts. Easy integration with tools for differentially methylated region (DMR) calling and annotation further eases the analysis of genome-wide methylation data. Overall, methrix enriches established WGBS workflows by bringing together computational efficiency and versatile functionality.

## Decoupling statistics and networks: a crowdsourced systematic assessment of transcription factor activity estimation from transcriptomics data

_Speaker_: Ricardo Omar Ramirez Flores

_Abstract_:
Transcriptome profiling followed by differential gene expression analysis often leads to lists of genes that are hard to analyze and interpret. Downstream analysis tools can be used to summarize deregulation events into a smaller set of biologically interpretable features. In particular, methods that estimate the activity of transcription factors (TFs) from gene expression are commonly used. It has been shown that the transcriptional targets of a TF yield a much more robust estimation of the TF activity than observing the expression of the TF itself. Consequently, for the estimation of transcription factor activities, a network of transcriptional regulation is required in combination with a statistical algorithm that summarizes the expression of the target genes into a single activity score. Over the years, many different regulatory networks and statistical algorithms have been developed, mostly in a fixed combination of one network and one algorithm. To systematically evaluate both networks and algorithms, we developed decoupleR (https://github.com/saezlab/decoupleR), an R package that allows users to apply efficiently any combination provided. We benchmarked an initial set of 84 combinations, comprising 7 methods and 13 networks. We evaluated the precision of different combinations in recovering perturbed TFs from different collections of gene expression datasets. Additionally, we tested the effects of combining multiple sources and estimations. We set up the package in a modular way which makes it easy and intuitive to extend it with further statistics or networks. We invite the community to participate by implementing their own statistics or integrating their gene regulatory network. With the decoupleR package, we lay the foundation for a crowdsourced systematic assessment of transcription factor activity estimation from transcriptomics data.

## DEWSeq: a new data analysis package for eCLIP data showing improved RNA-binding site assignments for the ENCODE eCLIP data

_Speaker_: Thomas Schwarzl

_Abstract_:
Recently, more than 3.000 human RNA-binding proteins (RBPs) have been discovered. Enhanced Crosslinking and Immunoprecipitation (eCLIP) sequencing has become a popular method for the detection of the RNA targets of RBPs. The ENCODE consortium has reported 223 eCLIP data sets and provided a bioinformatic analysis workflow for target site detection. To overcome limitations of the ENCODE analysis workflow regarding the specificity and reproducibility of assigned target sites, we developed a new analysis approach using the custom-built Python package “htseq-clip” for preprocessing, and our Bioconductor/R package DEWSeq for sliding-window significance testing with DESeq2 and binding site post-processing. We reanalyzed all 223 ENCODE eCLIP data sets and demonstrate a significant improvement in binding site detection based on a benchmark of known binding motifs. We performed functional analyses on these revised binding sites and disease associations from OpenTargets leading to new insights into disease networks, RBP regulatory feedback loops, and motifs.

## dasper: Detection of aberrant splicing events from RNA-sequencing data

_Speaker_: David Zhang

_Abstract_:
Although next-generation sequencing technologies have accelerated the discovery of novel gene-to-disease associations, the majority of Mendelian disease patients still leave the clinic without a genetic diagnosis. An estimated one third of these patients will have disorders caused by mutations impacting splicing. RNA-sequencing has been shown to be a promising diagnostic tool, however few methods have been developed to integrate RNA-sequencing into the diagnostic pipeline. Here, we introduce dasper, a Bioconductor package that improves upon existing tools for detecting aberrant splicing by using machine learning to incorporate disruptions in junction counts as well as coverage. dasper is designed for diagnostics, providing a rank-based report of how aberrant each splicing event looks, as well as including visualization functionality to facilitate interpretation (sashimi plots). We validate dasper using 16 fibroblast samples derived from patients with mutations known to impact splicing, and a further 50 patient samples as controls. We find that dasper is able to detect pathogenic splicing events with greater ease than existing leafcutterMD or z-score approaches. Furthermore, by only applying a broad OMIM gene filter (without any variant-level filters), dasper is able to detect pathogenic splicing events in the top 5 (median: 3.75) most aberrant splicing events for each patient. We also investigate the use of 504 GTEx fibroblast samples as controls and find that dasper leverages publicly available data effectively, ranking pathogenic splicing events in the top 25 (median: 24.25). Together, we believe dasper can increase diagnostic yield and enable the efficient implementation of RNA-sequencing for diagnostics in clinical laboratories.

## annoFuse and shinyFuse: Annotating, prioritizing, and interactively exploring putative oncogenic RNA fusions

_Speaker_: Federico Marini

_Abstract_:
Gene fusion events are a significant source of somatic variation across adult and pediatric cancers, and have provided some of the most effective clinically-relevant therapeutic targets. Still, computational algorithms for fusion detection from RNA sequencing data show low overlap of predictions across methods. Additionally, events such as polymerase read-throughs, mis-mapping due to gene homology, and fusions occurring in healthy normal tissue require stringent filtering, making it difficult for researchers and clinicians to discern gene fusions that might be true underlying oncogenic drivers of a tumor and in some cases, appropriate targets for therapy. Here, we present annoFuse and shinyFuse, an R package with a companion Shiny application, developed to annotate, identify, and explore biologically and clinically-relevant expressed gene fusions in single samples, as well as recurrent fusions in a given cohort. We applied annoFuse to STAR-Fusion and Arriba results for 1,028 pediatric brain tumor samples provided as part of the Open Pediatric Brain Tumor Atlas (OpenPBTA) Project. First, we used FusionAnnotator to identify and filter “red flag” fusions found in healthy tissues or in gene homology databases. Using annoFuse, we filtered out fusions known to be artifactual and retained high-quality fusion calls. Second, we prioritized and captured known and putative oncogenic driver fusions previously reported in TCGA or fusions containing gene partners that are known oncogenes, tumor suppressor genes, genes documented in the COSMIC database, or transcription factors. Finally, we enable seamless exploration of the output from the annoFuse workflow by providing a web application, shinyFuse, where users can explore fusion calls in interactive tables augmented with links to external databases, visualize breakpoints by transcript with domain annotation, and plot recurrent fusions within cohorts. Fusions are often associated with early tumorigenic events, therefore we anticipate that systematically prioritizing fusion data with annoFuse will help inform experimental and therapeutic studies. The annoFuse package is available at https://github.com/d3b-center/annoFuse, and we expect to submit it to the Bioconductor project.

## bambu: reference-guided transcript discovery and quantification for long read RNA-Seq data

_Speaker_: Ying Chen

_Abstract_:
Understanding the transcriptome is key in interpreting the function of the genome in human diseases. Quantification of transcriptome expression with short read sequencing, however, remains challenging because of isoform similarities and high repetitiveness across genomic regions. Long read Oxford Nanopore RNA sequencing can reduce the complexity of transcriptome quantification with reads that cover more than 80% length of the isoforms. However, this technology has a high sequencing error rate and often generates short fragmented reads due to RNA degradation or polyA homopolyphisms. Here, we present bambu, a long read isoform reconstruction and quantification method that performs probabilistic read assignment after extending annotations across samples to improve the accuracy of the transcriptome expression. We apply bambu to cancer cell line data with sequin spike-in controls and compare our method with estimates obtained from short read data. bambu successfully recovers unannotated isoforms,showed consistency in gene expression estimation with existing methods for short read RNA-Seq data, but improved accuracy in transcript expression estimation. Estimates from bambu can be used directly by other software (e.g., DESeq2) for downstream analysis to provide biological insights for studies in human diseases. bambu is implemented as a R package with an example data package NanoporeRNASeq, both available in Bioc 3.12.

## Inference of transcriptional and post-transcriptional dynamics from sequencing experiments

_Speaker_: Mattia Furlan

_Abstract_:
The combined action of RNA synthesis, processing and degradation, governed by the corresponding kinetic rates, determines genes expression levels and how they are modulated in response to a stimulus. The standard approach to study RNA metabolism at this resolution requires the profiling of both total and nascent transcription, and the inference of RNA kinetic rates based on a mathematical model of the RNA life cycle. INSPEcT is a Bioconductor package, available since 2015 (de Pretis S. et al, Bioinformatics 2015), which allows this analysis both at steady state and in time course. The profiling of nascent RNA is affected by a number of pitfalls. Therefore, we recently extended our tool with a novel approach (INSPEcT-) to estimate the kinetic rates from total RNA-seq data only (Furlan M. Et al, Genome Research 2020). This novel method provides consistent results and improves the inference of transcripts half-lives despite the reduction in experimental costs and complexity. It also permits the reanalysis of published data, and the study of RNA metabolism at unprecedented resolution for biological systems in which nascent RNA cannot be readily profiled. A little experience in mathematical modelling and coding is required to interpret and manipulate INSPEcT results, potentially reducing the number of users. For this reason, we recently developed a Shiny interface (INSPEcT-GUI) to enable researchers without specific expertise to easily interact with INSPEcT results and to explore how the kinetic rates shape gene expression programs (de Pretis S. et al, Frontiers in Genetics 2020). My contribution to the European Bioconductor Meeting concerns the overview of the whole INSPEcT analysis framework with a focus on INSPEcT-, its application to actual research questions, and INSPEcT-GUI.

## Correcting for Cell RNA Fractions in MDD RNA-seq Data

_Speaker_: Louise Huuki

_Abstract_:
Major Depressive Disorder (MDD) is a common condition characterized by periods of low mood, loss of interest, and low energy among other side effects. MDD is a heritable disease with complex genetic effects such as other conditions such as Bipolar Disorder. To further study these disorders, we will examine differences in transcriptional mechanisms between the brain tissue of individuals with MDD, Bipolar Disorder, and neurotypical control individuals. Our data consist of 1091 bulk RNA sequencing (RNA-seq) of postmortem human samples of the subgenual anterior cingulate cortex (sACC, n = 551) and the amygdala (n = 540) from all three diagnosis groups. However to ensure that downstream analysis of this RNA-seq data are only revealing differences between diagnosis, we must first control for RNA quality and differences in cell-type RNA fractions. Using the deconvolution algorithm MuSiC (Wang et al, Nat. Comms., 2019), with complementary snRNA seq data (Tran et al, bioRxiv, 2020) we have predicted the RNA fraction for six broad cell types for each sample. We have also observed significant correlation between cell-type RNA fractions and quality surrogate variables (Jaffe et al, PNAS, 2017) that are used for controlling RNA degradation effects. These results will inform what statistical model to use to explore differences among the diagnosis in the sACC and amygdala in the human brain.

## SourceSet: a graphical model approach to identify primary genes in perturbed biological pathways

_Speaker_: Elisa Salviato

_Abstract_:
 INTRODUCTION: The Gene Set Analysis (GSA) is one of the most used approaches for analyzing gene expression data: it moves from a gene-centered perspective towards a gene set-centered perspective, where the gene set are usually defined as groups of functionally related genes. Among GSA tools, Topological Pathway Analysis (TPA) aims to improve inferential analysis by exploiting the explicit biological network information. Although numerous TPA methods for identifying dysregulated genes have been proposed, they use marginal approaches and are therefore unable to distinguish between the so-called "primary genes" representing the source of perturbation – for example, a mutation, a copy number variation or an epigenetic change – from those that are merely affected by the propagation of that perturbation. METHODS: We presented a new method, called SourceSet able to distinguish between the primary and the secondary dysregulation within a perturbed pathway. The proposed method compares gene expression profiles in two conditions and detects the differences in both the mean and the covariance parameters. Set within the framework of Gaussian graphical models, it compares all marginal and conditional distributions induced by the underlying graph, and uses the results to infer the set of primary genes potentially responsible for the differential behaviors. To make it applicable when the number of samples/replicates is far below the number of genes – a scenario that characterizes omics data – we adopted an ad-hoc ridge strategy for estimating the covariance matrix and a permutation approach. RESULTS: The method is implemented in the SourceSet R package (available on CRAN), which contains statistics and graphical devices aiding the user in interpreting the obtained results. It consists of six core functions that, given a list of pathways to be analyzed and a gene expression matrix, allow the user to: i) identify a source set and a secondary set of each graph; ii) pool results from single-pathway analyses to gain a global view of results and obtain replicable summaries of research findings through additional visualization tools and statistics; iii) connect with Cytoscape software environment to visualize, explore and manipulate chosen pathways in a dynamic manner. CONCLUSIONS: Extensive simulations show that, when a dysregulation is present, SourceSet demonstrates high sensitivity and specificity in all considered scenarios, even with a low number of samples. Apart from array-based gene expression, we successfully applied SourceSet to log counts or RPKM/FPKM in next-generation sequencing experiments. Other possible applications include protein abundances and metabolomic data.

## Estimation of transcription factor and pathway activities from bulk and single-cell transcriptomics data with dorothea and progeny

_Speaker_: Christian Holland

_Abstract_:
Transcriptome profiling followed by differential gene expression analysis often leads to lists of genes that are hard to analyze and interpret. Functional genomics tools are powerful approaches for downstream analysis, as they summarize the large and noisy gene expression space into a smaller number of biologically meaningful features. In particular, methods that estimate the activity of transcription factors and pathways from gene expression are popular. Here we present the tools DoRothEA and PROGENy that allow inferring the activity of these molecule classes, respectively. Both rely on gene sets comprising downstream affected genes instead of mapping transcripts level to process members. While these tools were originally developed for application in human data it has been shown that they can functionally characterize mouse data as well. With the emergence of single-cell RNA-seq data, these tools were benchmarked for the application in single-cell. Both tools i) are robust against low gene coverage, ii) are able to detect perturbed transcription factors/pathways, iii) preserve cell type-specific information while reducing noise in parallel, and iv) are biologically meaningful. We made those tools easily accessible for the community as Bioconductor packages.

## Switching between mendelian genes space and clinical phenotypes space contributes to comparative Mendelian gene sets analysis

_Speaker_: Alejandro Cisterna García

_Abstract_:
Introduction Diseases with a genetic basis are usually diagnosed by looking for causal mutations in a panel of genes specifically associated with the disease. Gathering all phenotypes associated with the genes in a panel delivers a general phenotype-level description beyond the disease under study. For the purpose of improving genetic diagnosis, we need methods to evaluate the gene level phenotypic similarity between diseases and to identify differential phenotypes between the gene sets associated. Comparative phenotypic similarity analysis between gene panels can demonstrate common pathophysiology and help to create genetic links between diseases through their gene sets. PhenoExam works on this principle by integrating databases like The Human Phenotype Ontology (HPO) and The Mouse Genome Database (MGD). PhenoExamWeb is an R package that performs (1) phenotype enrichment analysis on a gene set, (2) measures statistically significant phenotype similarities between gene sets and (3) detects significant differential phenotypes for them. Phenotypic Similarity between two groups of genes is performed by assessing the statistical significance of the Phenotypic Overlap Ratio (POR) between those (i.e. the number of common phenotypes between the gene sets). PhenoExamWeb uses the HPO, MGD, and CRISPRbrain databases. Results We evaluated the specificity of PhenoExamWeb in detecting phenotype similarities between gene sets by comparing genetic forms of epilepsy (261 genes from NIMGenetics epilepsy panel) and artificial gene sets constructed with variable POR with the original epilepsy gene set and additional genes and with similar phenotypic connectivity not associated to epilepsy. We performed 1000 simulations for different proportions of epilepsy genes (80%, 75%, 60%, 40%). We calculated the POR significance test between the real and the artificial gene sets and we assessed that PhenoExamWeb can distinguish well amongst the simulated gene sets with very similar phenotypes (P < 10-5). We then applied PhenoExam to the detection of differential phenotypes between gene sets by comparing two genetic diseases with similar symptoms: juvenile Parkinson disease (PD) with 35 genes and early onset dystonia (EOD) with 50 genes from Genomics England PanelApp. They shared 127 significant phenotypic terms (out of 236 unique significant phenotypic terms in both), that yields a POR of 0.538 (P < 0.001). Phenotype relevance association analysis for PD and EOD (i.e. whether the shared phenotypes are similar in relevance, i.e. in the number of genes associated with them, within each gene set) results in an adjusted R squared of 0.664 (P < 1.22x10-31) which suggests that an important portion of the common phenotypes are similar in relevance. The p-values were obtained through randomization of 1000 random gene sets. We actually see they share human phenotypic terms such as Tremor (HP:0001337), Bradykinesia (HP:0002067), Rigidity (HP:0002063), Dystonia (HP:0001332), Abnormal gait (MP:0001406) and mammalian terms like Neuron degeneration (MP:0003224). But we also detect differential phenotypes. For example, significant only in PD phenotypes include Astrocytosis (MP:0003354; P < 4.69x10-12), Substantia nigra gliosis (HP:0011960; P < 4.13x10-11), Neuronal loss in central nervous system (HP:0002529; P < 3.71x10-6) and Orthostatic hypotension due to autonomic dysfunction (HP:0004926; P < 9.93x10-6). On the other hand, we found significant only phenotypes in EOD such as Writer's cramp (HP:0002356; P < 1.36x10-9), Hypoplasia of the corpus callosum (HP:0002079; P < 2.96x10-5), Acanthocytosis (HP:0001927; P < 2.76x10-3), Microcephaly (HP:0000252; P < 2.43x10-4), Intellectual disability, mild (HP:0001256; P < 4.27x10-3) and Hyperactive deep tendon reflexes (HP:0006801; P < 4.91x10-3). Discussion PhenoExamWeb allows us to switch from the gene space and the phenotype space. With PhenoExamWeb we can identify the statistically significant phenotypes of a gene set. It is useful to distinguish between gene sets or diseases with very similar phenotypes through projecting genes into their annotation based phenotypical spaces. PhenoExamWeb effectively discovers links between phenotypic terms across annotation databases by integrating different annotation databases. All these findings are supported with interactive plots (see tutorials at GitHub project) to foster the visualization and interpretation of findings. With the PD and EOD example above, we clearly see they hold phenotype-level similarities but also potentially interesting differential phenotypes. Resources Github: github.com/alexcis95/PhenoExamWeb We are setting up a website with shiny and PhenoExam. We will have it ready for the conference.


## MicrobiomeExperiment: a new class for microbiome data

_Speaker_: Felix GM Ernst

_Abstract_:
Dedicated class structures are fundamental for the R/Bioconductor ecosystem. In microbiome research, the phyloseq class has become a commonly accepted standard. However, recent developments in Bioconductor classes and growing sample sizes have opened up new opportunities and requirements to extend this popular format in order to address emerging research needs. Extended support is needed for instance for hierarchically structured data, linking with other data types, and for general performance optimization. We propose MicrobiomeExperiment as a novel class structure for microbiome data. This extends the TreeSummarizedExperiment and SingleCellExperiment classes, which provide thoroughly tested tools for hierarchical data, spanning both sample and feature spaces. The new class adds support for additional feature information that is specifically relevant for microbiome experiments; it allows a seamless incorporation of more detailed sequence information based on existing classes such as the DNAStringSet(List). The MicrobiomeExperiment class inherits support for sparse matrices and multiple assays while providing improvements in speed and memory compared to the currently available solutions for microbiome data. Importantly, it provides conversion methods from the widely used phyloseq class as well as other raw data types, thus enabling seamless conversion. By linking microbiome data more tightly to other already established and strongly supported Bioconductor classes, we hope to reduce overlapping development efforts, improve the interoperability of available tools, and simplify the addition and ensure the long-term sustainability of the new tools. MicrobiomeExperiment is aimed at solving shortcomings in the current microbiome R ecosystem and it has the potential to become a widely accepted standard for microbiome data in the R/Bioconductor ecosystem. The class structure and the associated package ecosystem are now under active development as a joint effort of multiple research teams, and feedback from the community will be highly valued.

## ROBIN (ROBustness In Network): for Comparison and Validation of communities

_Speaker_: Valeria Policastro

_Abstract_:
Clustering is an unsupervised machine learning procedure which aims to identify groups of unlabelled objects such as genes, cells, molecules, patients according to a similarity measure. Indeed, for large datasets, community detection is a very flexible and scalable approach for which many algorithms have already been developed. Nevertheless, there still is a lack of methodologies able to define the reliability of these clustering algorithms. To help to fill this gap, we created robin (ROBustness In Network), an R/CRAN package based on a statistical method able to assess the robustness of the community structure found by one or more methods and to give indications about their reliability. The procedure firstly detects if the community structure found by a set of algorithms is statistically significant and secondly compares two selected detection algorithms on the same graph to choose the one that better fits the network of interest. Despite our approach can be used on any kind of graph-based clustering, here we show its applicability to clustering analysis on a single cell dataset.

## Assessment of statistical methods from single cell, bulk RNA-seq, and metagenomics applied to microbiome data

_Speaker_: Matteo Calgaro

_Abstract_:
The correct identification of differentially abundant microbial taxa between experimental conditions is a methodological and computational challenge. Recent work has produced methods to deal with the high sparsity and compositionality characteristic of microbiome data, but independent benchmarks comparing these to alternatives developed for RNA-seq data analysis are lacking. We compare methods developed for single-cell and bulk RNA-seq, and specifically for microbiome data, in terms of suitability of distributional assumptions, ability to control false discoveries, concordance, power, and correct identification of differentially abundant genera. We benchmark these methods using 100 manually curated datasets from 16S and whole metagenome shotgun sequencing. The multivariate and compositional methods developed specifically for microbiome analysis did not outperform univariate methods developed for differential expression analysis of RNA-seq data. We recommend a careful exploratory data analysis prior to application of any inferential model and we present a framework to help scientists make an informed choice of analysis methods in a dataset-specific manner.

## BNPmix: an R package to estimate Bayesian nonparametric mixtures

_Speaker_: Riccardo Corradin

_Abstract_:
Bayesian nonparametric mixtures are flexible models for density estimation and clustering analysis, composed by an infinite number of components. Different strategies were proposed in the last decades to deal with the estimation of these models, most resorting to Markov chain Monte Carlo (MCMC) methods. Along with the study of a new approach to estimate infinite mixtures, we implemented an efficient R package, named BNPmix. The package include three different MCMC strategies to deal with nonparametric mixtures, and allows for flexible estimation of the models by choosing the algorithm, and by tuning model/algorithm-specific parameters. Routines for post-processing the results are also available. In order to face the computational complexity of the problem, the main functions are written in C++ and interfaced with R through the Rcpp and RcppArmadillo packages. The BNPmix packages is integrated in the R environment for Bayesian analysis and graphical visualization.

## Browsing and searching the Bioconductor codebase

_Speaker_: Mike Smith

_Abstract_:
The Bioconductor project uses a central Git server (git.bioconductor.org) as the foundation of the package build system, and is the source for all code that is distributed to users. Read access to all packages is open universally, meaning anyone can clone a package and examine its source code. However there is no facility for browsing the Bioconductor code base in a graphic interface like those provided by GitHub or GitLab. This can be particularly problematic if one wants to search across all packages e.g. to find examples of how other developers have implemented something, or whether a particular function in one package is used by any other. BioC Code Tools addresses these issues by using existing open source tools. Firstly, it provides a graphical interface to all packages using P3X Gitlist (https://github.com/patrikx3/gitlist), giving easy access both to the current status of the package as well as to branches, commit history etc. Secondly, cross-package code searching is implemented using zoekt (https://github.com/google/zoekt). This flexible search engine allows the use of regular expressions to search file content as well filtering based file name or type, providing a mechanism to iteratively drill down through the entire codebase to find something.

## The Bioconductor teaching committee: A collaborative effort to consolidate Bioconductor-focused training material and establish a community of trainers

_Speaker_: Charlotte Soneson

_Abstract_:
The Carpentries (https://carpentries.org) aims to teach foundational computational and data science skills to researchers using local instructors that have been specially trained in sound pedagogical methods. They offer a variety of introductory R modules, both on generic scripting and geared towards particular fields like ecology and genomics, but there are none currently involving Bioconductor-focused training material. The Bioconductor teaching committee was established in early 2020, with the aims to: (1) provide networking opportunities and coordinate training activities (e.g., Carpentries instructor certification) for members of the Bioconductor community interested in education, (2) establish a connection with the Carpentries and begin assembling introductory as well as more specialized Bioconductor-focused training material consistent with Carpentries guiding principles, and (3) coordinate delivery of the developed training material at Bioconductor-related events, such as conferences or standalone courses, both in person and virtually. The group meets regularly and is open to interested members of the Bioconductor community. In this short talk, we will provide an overview of the current state of our efforts, and invite a more general discussion around the training needs and opportunities within the Bioconductor community and potential for improvements through collaboration with The Carpentries.

## "Empower" the Biologist and "liberate" the Bioinformatician

_Speaker_: Naji Faris

_Abstract_:
This talk and demo are about how to “empower” the Biologist and “liberate” the bioinformatician for data analysis using the platform called Tercen. Tercen supplies a powerful visual and relational layer to the Bioconductor R packages, allowing Bioconductor packages to be accessed by non-coders. Tercen aims to make biologists relax with data, even if there is a deluge. Tercen helps non-coders get back control of their research using bioconductor packages for clinical data, flow cyto, masscyto, RNAseq, single cell RNAseq, mass spec, plate readers, and more. Tercen allows biologists to link a sequence of Bioconductor packages in innovative and synergistic ways. Tercen aims to liberate bioinformaticians by removing them from the operational support and allow them to concentrate on the algorithmic fundamentals of Bioconductor R packages. Tercen is an early stage European startup in the area of open science. The talk covers use cases where biologists and bioinformaticians have benefited from the approach. Feedback from the audience is very much appreciated and will influence the roadmap of the Tercen platform.

## Interactive and comprehensive exploration of ExpressionSet using ExpressionSetViewer

_Speaker_: Chen Meng

_Abstract_:
ExpressionSet is a standard S4 object storing high-throughput omics data in Bioconductor. It consists mainly of three components 1) the expression matrix where the rows are features (genes, mRNAs, proteins, etc) and columns are the samples (cell lines, patients, etc); 2) phenotype data - the annotation information of samples and 3) feature data – the annotation information of features. Interpreting biology from omics data relies on the integrative analysis of the data triplet and a wide range of statistical methods are often used to answer the different biological questions. Here, I describe a new package “ExpressionSetViewer”, which visualizes ExpressionSet in an interactive way. The ExpressionSetViewer has a separate back- and front-end. In the back-end, users need to prepare an ExpressionSet that contains all the necessary information for the downstream data interpretation. Therefore, we imposed some extra requirement on the headers of phenotype data or feature data so that the provided information can be clearly recognized by the front-end, at the same time, keep a minimum modification on the existing ExpressionSet object. The pure dependency on R/Bioconductor guarantees maximum flexibility in the statistical analysis in the back-end. Once the expression set is prepared, it can be visualized using the front-end, implemented by shiny and plotly. Both features and samples could be selected from (data) tables or graphs (scatter plot/heatmap). Different types of analyses, such as enrichment analysis (using Bioconductor package fgsea or fisher's exact test) and STRING network analysis, will be performed on-the-fly and the results are visualized simultaneously. When a subset of samples and a phenotype variable is selected, a significance test on means (t-test or ranked based test; when phenotype variable is quantitative) or test of independence (chi-square or fisher’s exact test; when phenotype data is categorical) will be performed to test the association between the phenotype of interest with the selected samples. Therefore, ExpressionSetViewer will greatly facilitate data exploration, many different hypotheses can be explored in a short time without the need for knowledge of R. In addition, the resulted data could be easily shared using a shiny server. Otherwise, a standalone version of ExpressionSetViewer together with designated omics data could be easily created by integrating it with portable R, which can be shared with collaborators or submitted as supplementary data together with a manuscript.

## Scalable or density-preserving (but not both) non-linear dimensionality reduction with snifter and densvis

_Speaker_: Alan Brian O'Callaghan

_Abstract_:
The dimensionality reduction techniques t-distributed stochastic neighbor embedding (t-SNE) and uniform manifold approximation and projection (UMAP) are widely used for visualizing single-cell RNA-sequencing (scRNA-seq) data. These algorithms allow for the representation of complex high-dimensional transcriptomic data in 2-dimensions. This is often useful for interpreting the results of complex analyses on large datasets, such as the identification of cell lineages or distinct cellular subtypes or subgroups. However, t-SNE scales poorly to large datasets, and both t-SNE and UMAP can produce misleading plots by failing to accurately represent the varying transcriptional heterogeneity within the original high-dimensional space. We present densvis and snifter, two R packages that attempt to deal with these issues to improve the quality and scalability of these dimensionality reduction techniques in the Bioconductor ecosystem.

## Peptide Correlation Analysis (PeCorA) Reveals Differential Proteoform Regulation

_Speaker_: Maria Dermit

_Abstract_:
Shotgun proteomics techniques infer the presence and quantity of proteins using peptide proxies produced by cleavage of the proteome by a protease. Most protein quantitation strategies assume that multiple peptides derived from a protein will behave quantitatively similar across treatment groups, but this assumption may be false due to (1) heterogeneous proteoforms and (2) technical artifacts. Here, we describe a strategy called peptide correlation analysis (PeCorA) that detects quantitative disagreements between peptides mapped to the same protein. PeCorA fits linear models to assess whether a peptide’s change across treatment groups differs from all other peptides assigned to the same protein. PeCorA revealed that ~15% of proteins contain at least one discordant peptide. Inspection of the discordant peptides shows utility of PeCorA for direct and indirect detection of regulated PTMs, and also for discovery of poorly quantified peptides. Exclusion of poorly quantified peptides before protein quantity summarization decreased false positives in a benchmark dataset. Finally, PeCorA of human plasma proteomics revealed lower abundance of one inactive isoform of coagulation enzyme thrombin specific to COVID-19 patients. PeCorA is freely available as an R package that works with arbitrary tables of quantified peptides.

## Prolfqua - R package for proteomics label-free quantification

_Speaker_: Witold Wolski

_Abstract_:
The package for proteomics label free quantification prolfqua (read: prolevka) evolved from functions and code snippets used to visualize and analyze label free quantification data. To compute protein fold changes among treatment conditions, we first used t-test or linear models and then used functions implemented in the package limma to obtain moderated p-values. We evaluated MSStats, ROPECA or MSqRob all implemented in R, with the idea to integrate the various approaches. Although all these packages are written in R, model specification, input and output formats among them differ widely and wildly, which made our first attempt to use the original implementations challenging. Therefore, and also to understand the algorithms used, we attempted to reimplement those methods, where possible. The R-package prolfqua is the outcome of this venture. When developing prolfqua, we draw inspiration from packages such as sf, which use data in a long table format, dplyr for data transformation and ggplot2 for visualization. In the long table format each column stores a different attribute, e.g.~there is only a single column with the intensities. In prolfqua the data needed for analysis is represented using a single data-frame in long format and an R6 configuration object. The configuration annotates the table, i.e.~specifies what information is in which column. The use of an annotated table makes integrating new data if provided in long formatted tables simple. Therefore, all that is needed to incorporate Spectronaut, Skyline text output, or MSStats inputs is to update the configuration object. For software like MaxQuant writing the data in a wide table format, with several intensity columns, one for each sample, we implemented methods that transform the data into a long format. Relying on the long data table format enabled us to easily access a large variety of useful data manipulation and visualizations methods implemented in the R packages dplyr and ggplot2. A further design decision, which sets prolfqua apart is that it embraces R's linear model formula interface, including the lme4 mixed effect models formula interface. R's formula interface for linear models is flexible, widely used and well documented. These interfaces allow specifying a wide range of essential models including parallel designs, factorial designs, repeated measurements and many more. Since prolfqua uses R modelling infrastructure directly, we can fit all these models to proteomics data. This is not easily possible with any other package dedicated to proteomics data analysis. For instance, MSStats, although using the same modelling infrastructure, supports only a subset of possible models. Limma, supports the R formula interface but not for linear mixed models. Since ROPECA relies on limma it is limited to the same set of models. MSqRob allows specifying fixed and random effects, however, it is unclear how interactions among factors can be specified, estimated or tested. The use of R's formula interface does not limit prolfqua to the output provided by the R modelling infrastructure. prolefqua implements p-value moderations, and computes probabilities of differential regulation, as suggested by ROPECA. Last but not least, ANOVA analysis or model selection using the likelihood ratio test for thousand of proteins can also be performed. To use prolfqua knowledge of the R regression model infrastructure is of advantage. Acknowledging the complexity of the formula interface, we provide an MSstats emulator, which derives the model formula from an MSstats formatted input file. We benchmarked all the methods implemented in prolfqua: linear models, mixed effect models, p-value moderation, ROPECA, as well as Bayesian regression models implemented in brms using a benchmark dataset, enabling us to evaluate the practical relevance of these methods. Last but not least, prolfqua supports elements of the LFQ data analysis workflow, e.g.~computing coefficients of Variations (CV) for peptide and proteins, sample size estimation, visualization and summarization of missing data, summarization of intensities, multivariate analysis, etc. It also implements various protein intensity summarization and inference methods, e.g.~top 3, or Tukeys median polish. Our package makes it relatively easy to perform proteomics data analysis. We continue extening the functionality of the package, i.e.~currently we are adding methods for count data modelling. The package can be installed from www.github.com/wolski/prolfqua.

## CyTOF data analysis workflow with Bioconductor packages

_Speaker_: Anne-Maud Ferreira

_Abstract_:
Mass cytometry or Cytometry by Time Of Flight (CyTOF) is used to study cell population variation between different experimental conditions. This technology measures the expression of more than 40 proteins at single cell resolution by using heavy metal isotope labelled antibodies. Although the measurement technology is well established, CyTOF data analytic workflow methods are still expanding. Bioconductor provides a panel of methods and tools designed to analyze this high-throughput data. Using a selection of Bioconductor packages, we have studied the impact of the dengue virus (DENV) on natural killer (NK) cells during pediatric and adult infections. The FlowCore package allows importing the data in R. The SingleCellExperiment package provides a core data structure and framework for downstream analysis. Integrated into the CATALYST package, the combination of the flowSOM and ConsensusClusterPlus algorithms allow clustering the cells for cell phenotyping. Using this workflow, we identified different subsets within the NK cell population. With the help of the diffcyt package, we performed differential abundance tests between healthy and DENV positive patients within the pediatric and adult populations. We found that pediatric DENV infection leads to changes in the total NK cell frequency within two NK cells subsets. Next, we used a custom-made package CytoGLMM to identify different proteins within adult and pediatric samples that are predictive of the DENV infection. The results of these analyses suggest that the disease severity in pediatric dengue cases may be partially explained by a diminution in the NK cell activation.

## POMA: An user-friendly workflow for pre-processing and statistical analysis of mass spectrometry data

_Speaker_: Pol Castellano Escuder

_Abstract_:
Mass spectrometry, like other high-throughput technologies, usually faces a data mining challenge to provide an understandable output to advance in biomarker discovery and precision medicine. Often, statistical analysis is one of the hard points and it’s critical in the subsequent biological interpretation of the results. Due to this fact combined with the computational programming skills needed for this type of analysis, several bioinformatic tools have emerged to simplify mass spectrometry data analysis. However, sometimes the analysis is still limited to a few hidebound statistical methods and to a low-flexible datasets. POMA introduces a structured, reproducible and easy-to-use workflow for the visualization, pre-processing, exploratory and statistical analysis of mass spectrometry data. In summary, POMA enables a flexible data cleaning and statistical analysis in one comprehensible and user-friendly R/Bioconductor package (https://bioconductor.org/packages/release/bioc/html/POMA.html). POMA also has a Shiny app version with all the package functions implemented. See https://github.com/pcastellanoescuder/POMAShiny/.

# Posters

## Interactions between explicit and latent covariates detect cell-type specific treatment effects in exploratory single cell analysis without clustering

_Presenter_: Constantin Ahlmann-Eltze

_Abstract_:
Single cell studies increasingly move to replicated experimental designs to investigate cell type specific effects and not just document tissue heterogeneity. However, existing tools for exploratory data analysis struggle to incorporate known covariates into the analysis. Here, we present a new framework to incorporate known covariates and identify non-linear effects. We extend the classical the classical linear model / factor analysis framework with interaction terms between explicit and latent covariates. We identify interpretable factors that describe how subpopulations react differentially to treatment and are able to correct for known batch effects.

## Bioconductor packages for analysis of thermal proteome profiling data

_Presenter_: Nils Kurzawa

_Abstract_:
Thermal proteome profiling (TPP) is a mass-spectrometry-based technique originally developed for detection of drug targets through compound-binding induced shifts in thermal stability of proteins. In the last years, however, it has been realized that TPP can also be applied to study low-affinity protein interaction partners such as metabolites and is also able to inform on protein-protein interactions. Bioconductor packages for analyzing different TPP assay formats and aspects of the data have been implemented in parallel to the continuous development of the experimental assays. Here, I would like to introduce available Bioconductor packages for analysis of TPP experiments and briefly outline their functionality.

## HCAMatrixBrowser: Making use of the HCA Matrix API in Bioconductor

_Presenter_: Marcel Ramos

_Abstract_:
The Human Cell Atlas Data Coordination Platform (DCP) provides queryable matrices via a REST API. The $HCAMatrixBrowser$ makes use of the rapiclient to send queries to the HCA matrix endpoint and download expression data. The package provides standard Bioconductor representations based on the compatible formats from the API. It uses the LoomExperiment package to serve matrix data from an HDF5 loom format and supports the matrix market format (MTX) conversion to SingleCellExperiment.

## Activation of digestive-tract and stratified-epithelium transcriptional programs by promoter DNA demethylation in human tumor cells

_Presenter_: Anna Diacofotakis

_Abstract_:
Our lab has previously shown that a wide variety of human tumors exhibit aberrant activation of genes that are normally expressed specifically in testicular germ cells. Their activation in tumors is due to the demethylation of their promoter region, a consequence of genome-wide DNA methylation loss in these cells. It is now recognized that several of these ‘cancer-germline’ genes contribute to tumor development. To date, DNA demethylation in tumors has essentially been associated with the activation of germline-specific genes. Whether this epigenetic alteration induces other tissue-specific programs in tumors remains an open question. By combining RNA-seq data and methylation data of a series of lung adenocarcinoma cell lines (LUAD), we now show that two additional tissue-specific expression programs are aberrantly activated in tumors due to promoter DNA demethylation. The first group of genes was found to be physiologically expressed in different gastro-intestinal tissues such as stomach, liver, small intestine or colon. We termed this group of genes ‘lower digestive tract genes’ (LDG). The second group of genes was more intriguing as its expression was restrained to three seemingly unrelated tissues: esophagus, skin and vagina. However, a common feature of these tissues is that they harbour a stratified squamous epithelium. We therefore refer to this group of genes as ‘stratified squamous epithelial’ genes (SSE). Using TCGA project, we confirmed that these two transcriptional programs are activated in LUAD tumors and expanded this observation to a variety of epithelial tumors. We verified experimentally the dependence of these LDG and SSE programs on DNA methylation for their regulation by treating a series of cell lines with DNA methylation inhibiting agent 5-azadeoxycytidine. Furthermore, we found that to the contrary of cancer-germline genes, LDG and SSE genes are activated in tumors irrespective of the overall DNA methylation content of these tumors. This indicates that these two group of genes are not activated as a consequence of global genome hypomethylation, but rather through targeted promoter demethylation processes. Finally, we show that the activation of some genes belonging to the SSE transcriptional program is associated with clinical parameters in LUAD such as a higher tumor grade and poor patient survival. Overall, our study identified two transcriptional circuits that are dependent on DNA demethylation for their activation in tumors. Our data point to an activation mechanism involving targeted DNA demethylation events in these cells. Finally, the activation of these genes could be potentially exploited in the clinic as prognostic markers of LUAD.

## CoExp Web Application, a web application for the interactive use of co-expression networks

_Presenter_: Sonia García Ruiz

_Abstract_:
Gene co-expression network analysis (GCNA) is a powerful tool to reveal linkages between clusters of highly co-expressed genes and their potentially unknown biological functions. Several R software packages, such as the km2gcn or WGCNA, have provided a successful mean to uncover these underlying relationships. However, whereas the interest in exploring these genetic interconnections may arise from people presenting diverse backgrounds and levels of expertise, the path leading to such exploitation generally requires a minimum degree of coding skills. In addition, the general tendency to represent co-expression networks through undirected graphs makes the development of a web-based platform for the implementation of GCNA highly appropriate. Here we present CoExp Web Application, a user-friendly web-based tool that offers the exploitation of 109 co-expression networks made available through the CoExpNets R package. CoExp Web Application allows (1) exploring the different networks by the ontology terms and cell types that more strongly correlate with their gene clusters; (2) annotating a list of preferred genes within the modules that belong to a collection of selected networks; (3) generating an undirected graph to represent the most important genes neighbouring a gene of interest; and (4) downloading the raw data used in all the beforementioned operations for third-party usage. Finally, CoExp Web Application has been configured to be used from a wide range of devices enabling web browsing as well as published through the Docker container technology, offering the possibility of a local installation and use of the tool.

## SingleCellMultiModal: Integrating Multi-modal Single Cell Experiment datasets

_Presenter_: Kelly Eckenrode

_Abstract_:
Developments in sequencing and computational technology have created a surge of single cell multimodal experimental methods which assay complementary molecular data on the same individual cells. However, analysis of these data is hindered by complex and non-standard data representation, manipulation, and analysis. We present a new Bioconductor package, SingleCellMultiModal, that provides ready-to-use datasets from mouse and human samples for benchmarking and methods development from single-cell multimodal technologies including CITE-seq, ECCITE-seq, scNMT, M&T, and seqFISH. SingleCellMultiModal provides pre-processed datasets using the integrative MultiAssayExperiment Bioconductor class to coordinate management of multiple experimental assays performed on an overlapping set of cells and samples. The package is intended to facilitate development of bioinformatic and statistical methods in Bioconductor to meet the challenges of tracking molecular layers to phenotypic outputs including cell differentiation, activity, and disease.

## GenomicSuperSignature: interpretation of RNA-seq experiments through robust, efficient comparison to public databases

_Presenter_: Sehyun Oh

_Abstract_:
PURPOSE: Thousands of RNA sequencing profiles have been deposited in public archives, yet remain unused for the interpretation of most newly performed experiments. Methods for leveraging these public resources have focused on the interpretation of existing data, or analysis of new datasets independently, but do not facilitate direct comparison of new to existing experiments. The interpretability of common unsupervised analysis methods such as Principal Component Analysis would be enhanced by efficient comparison of the results to previously published datasets. METHODS: To help identify replicable and interpretable axes of variation in any given gene expression dataset, we performed principal component analysis (PCA) on 536 studies comprising 44,890 RNA sequencing profiles. Sufficiently similar loading vectors, when compared across studies, were combined through simple averaging. We annotated the collection of resulting average loading vectors, which we call Replicable axes of Variation (RV), with details from the originating studies and gene set enrichment analysis. Functions to match PCA of new datasets to RVs from existing studies, extract interpretable annotations, and provide intuitive visualization, are implemented as the GenomicSuperSignature R package, to be submitted to Bioconductor. RESULTS: The process to generate RVs is robust to the batch effects and the presence of low-quality or irrelevant studies. We demonstrate that the GenomicSuperSignature package allows nearly instantaneous matching of PCA axes in new datasets to pre-computed RVs from hundreds of previous studies and thousands of samples on an ordinary laptop. GenomicSuperSignature thereby links, in negligible compute time, relevant studies from the public literature, and provides interpretation of potentially subtle biological signals through gene set enrichment analysis of loadings and MeSH term analysis of matching studies. We demonstrate RVs associated with phenotype can provide insight into weak or indirectly measured biological attributes in a new study by leveraging accumulated data from published datasets. The performance of GenomicSuperSignature is compared to complementary previous works, and demonstrates 1) ability to identify replicable PCA loadings for colorectal carcinoma using a large, non-specific RNA-seq database instead of a focused microarray database, and 2) estimate neutrophil counts through transfer learning on new data comparably to the previous efforts despite major differences in training datasets and model building processes. CONCLUSION: GenomicSuperSignature enables researchers to analyze new data in the context of existing databases with minimal computing resources.
