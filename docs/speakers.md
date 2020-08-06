---
layout: default
speakers:
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
---

{% include header.md %}

## Confirmed speakers

{% for s in page.speakers %}
{% assign imgpath = "images/speakers/" | append: s.name | remove: ' ' | append: '.jpg' %}
<img src="{{ imgpath }}" style="float:right; width:120px; height:150px; object-fit: cover">
### [{{ s.name }}]({{ s.url }}), {{ s.inst }}

> {{ s.blurb }}

{% endfor %}
