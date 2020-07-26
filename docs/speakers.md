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
---

{% include header.md %}

## Confirmed speakers

{% for s in page.speakers %}
{% assign imgpath = "images/speakers/" | append: s.name | remove: ' ' | append: '.jpg' %}
<img src="{{ imgpath }}" style="float:right; width:120px; height:150px; object-fit: cover">
### [{{ s.name }}]({{ s.url }}), {{ s.inst }}

> {{ s.blurb }}

{% endfor %}
