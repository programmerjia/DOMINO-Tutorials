# Welcome to DOMINO’s documentation!

# Overview
DOMINO is built based on a self-supervised multi-view graph contrastive learning framework. It is designed to integrate spatial coordinates and gene expression information for robust identification of tissue domains from spatial transcriptomics data. DOMINO employs graph neural networks (GNNs) as base encoder, constructing a multi-view graph contrastive learning framework using the original graph and the diffusion graph to learn spot representations in the ST data. After representation learning, spatial clustering is performed using clustering tools such as mclust to assign spots to corresponding spatial domain.


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

**The overall framework of DOMINO**
**a**. DOMINO uses spatial coordinates and gene expression of each cell from spatial transcriptomics data as model inputs to construct an undirected graph. A diffusion graph is further derived using graph diffusion convolution with personalized PageRank.
**b**. Multiple graph views are encoded through graph convolutional networks. The original graph, shuffled graph, and diffusion graph each generate node-level and graph-level embeddings. These embeddings are used to build multi-view representations in latent space.
**c**. A self-supervised contrastive learning framework is applied by defining positive and negative embedding pairs across different graph views. 
**d**. The learned low-dimensional embeddings can then be used to identify spatial domains, which can further be used in different downstream analyses, including cell type composition analysis, differential expression, and inference of cell–cell communication. 

