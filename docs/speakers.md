---
layout: default
speakers:
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
