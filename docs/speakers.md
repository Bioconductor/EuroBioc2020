---
layout: default
speakers:
    - name: Sehyun Oh
      inst: The City University of New York
      url: https://cunyisph.org/team/sehyun-oh/
      blurb: "Sehyun Oh is a postdoctoral fellow at the CUNY Graduate School of Public Health and Health Policy in Levi Waldron’s group. She is a molecular biologist by training with a PhD from University of Minnesota and a postdoc from Columbia University and transitioned into dry-lab bioinformatics. Her current project is focused on the development of approaches for multi-omic data analysis and implementation of bioinformatics workflows in the cloud-based platform."
    - name: Holger Heyn
      inst: CNAG-CRG
      url: https://www.cnag.crg.eu/teams/genome-research-unit/single-cell-genomics-team
      blurb: "Holger Heyn is the Team Leader of the Single Cell Genomics Group at the National Centre for Genomic Analysis (CNAG-CRG) in Barcelona. His team's mission is the implementation of latest single-cell and spatial sequencing technologies and their application in a basic research (Human Cell Atlas) and translational (Immuno-oncology) context."
    - name: Henrik Bengtsson
      inst: UCSF
      url: https://github.com/HenrikBengtsson
      blurb: "Henrik Bengtsson is the author of over 30 R packages in CRAN and Bioconductor, including matrixStats. His research is on statistics and bioinformatics with an emphasis on high-quality, reproducible method development, sustainable implementations, and large-scale processing. He is an Associate Professor in the Department of Epidemiology & Biostatistics at the University of California, San Francisco (UCSF), affiliated with the UCSF Helen Diller Family Comprehensive Cancer Center, and a member of the R Foundation and the R Consortium Infrastructure Steering Committee."
    - name: Monica Chiogna
      inst: University of Bologna
      url: https://www.unibo.it/sitoweb/monica.chiogna2
      blurb: "Monica Chiogna is Professor of Statistics at the University of Bologna. Her research interests are in the application of statistical methods to biological and medical sciences. Recent projects include evaluation of tests and biomarkers for disease screening, diagnosis, prognosis and risk prediction and the use of  graphical models as comprehensive probabilistic tools to model biological networks and to study their perturbations."
    - name: Luke Zappia
      inst: Helmholtz Zentrum München
      url: https://lazappi.id.au/
      twitter: https://twitter.com/_lazappi_
      github: https://github.com/lazappi
      blurb: "Luke Zappia is a postdoctoral researcher in the Theis lab at the Helmholtz Zentrum München Institute of Computational Biology and the Technische Universität München. His work focuses on methods for analysing scRNA-seq data and during his PhD with Alicia Oshlack he developed the Splatter Bioconductor package for simulating scRNA-seq data."
    - name: Britta Velten
      inst: German Cancer Research Center
      url: https://bv2.github.io/
      twitter: https://twitter.com/brittavelten
      github: https://github.com/bv2
      blurb: "Britta Velten is a postdoctoral researcher at the German Cancer Research Center. Originally trained as mathematician she gained her PhD from ETH Zurich and EMBL working with Peter Bühlmann and Wolfgang Huber. In her research, she aims to apply statistical reasoning in combination with modern machine learning approaches to gain insights into fundamental processes that underpin biological systems. For this she develops statistical methods and computational tools to analyse and integrate multi-omics data, with applications ranging from personalised medicine to single cell biology and developmental biology."
    - name: Elsa Bernard
      inst: Memorial Sloan Kettering Cancer Center
      url: https://elsab.github.io/
      blurb: "Elsa Bernard is a research fellow in Computational Oncology at Memorial Sloan Kettering Cancer Center in Elli Papaemmanuil's group. She obtained her PhD in Bioinformatics from the Center for Computational Biology at MinesParistech/Institut Curie under the supervision of Jean-Philippe Vert. Her work focuses on the development and application of statistics and machine learning techniques to cancer genomics with a focus on hematologic malignancies. Her current project includes the molecular characterization of >3000 patients with myelodysplastic syndromes and the development of personalized prognostic and predictive risk models."
    - name: Nicola Segata
      inst: University of Trento
      url: http://segatalab.cibio.unitn.it/
      blurb: "Nicola Segata, Ph.D., is Associate Professor at the CIBIO Department at the University of Trento (Italy). His lab comprises more than 20 researchers and employs experimental metagenomic tools and novel computational approaches to study the diversity of the microbiome across conditions and populations and its role in human diseases and infections. The projects in the lab bring together computer scientists, microbiologists, statisticians, and clinicians and are generally focused on profiling microbiomes with strain-level resolution and on the meta-analysis of very large sets of metagenomes with novel computational tools. His work is supported by the European Research Council and by several other European agencies."
---

{% include header.md %}

## Confirmed speakers

{% for s in page.speakers %}
{% assign imgpath = "images/speakers/" | append: s.name | remove: ' ' | append: '.jpg' %}
<img src="{{ imgpath }}" style="float:right; width:120px; height:150px; object-fit: cover">
### [{{ s.name }}]({{ s.url }}), {{ s.inst }}

> {{ s.blurb }}

{% endfor %}
