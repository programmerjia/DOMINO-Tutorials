# Welcome to DOMINO’s documentation!

# Identifying distinct spatial domains from clear cell and endometrioid ovarian carcinoma using DOMINO


```{toctree}
:maxdepth: 1
:caption: Documentation Contents
Installation
Tutorial 1
Tutorial 2
```

```{image} /_static/overview.png
:width: 800px
:alt: DOMINO Architecture Overview
:align: center
:class: shadow
```

# Overview
DOMINO is built based on a self-supervised multi-view graph contrastive learning framework. It is designed to integrate spatial coordinates and gene expression information for robust identification of tissue domains from spatial transcriptomics data. DOMINO employs graph neural networks (GNNs) as base encoder, constructing a multi-view graph contrastive learning framework using the original graph and the diffusion graph to learn spot representations in the ST data. After representation learning, spatial clustering is performed using clustering tools such as mclust to assign spots to corresponding spatial domain.