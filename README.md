📖 Project Overview
This project aims to analyze the spatial-level differences in gene expression within brain tissue. By leveraging Spatial Transcriptomics (Method of the Year 2020) and the GASTON deep learning method, we successfully mapped the topography of tissue slices and identified distinct functional layers through differential expression analysis.


🧬 Background
The brain exhibits significant spatial heterogeneity in gene expression. Traditional bulk sequencing loses this spatial context. We utilized the methodology described in the paper "Mapping the topography of spatial gene expression with interpretable deep learning" to reproduce the GASTON method. This approach allows for the identification of continuous spatial domains (isodepths) and gene expression gradients.


📂 Dataset
Source: Gene Expression Omnibus (GEO) 
Accession Number: GSE184510 
Specific Sample Used: GSM5591748 (HE_ON_2339) 


🛠 Methodology
Data Preprocessing: Input spatial coordinates and gene expression vectors.
GASTON Implementation: We reproduced the GASTON method to learn the "topography" of the tissue slice from spatial transcriptomics data.
Clustering: Post-GASTON processing involved clustering to unveil distinct brain layers (Layer 0, 1, 2, 3).
Differential Expression Analysis: We performed differential gene expression analysis across the identified layers to determine their biological functions and cell-type dominance.


📊 Key Results: Layer Characterization
Our analysis identified four distinct layers with specific gene expression signatures:

Layer 3: White Matter / Oligodendrocytes
Characteristics: Dominated by oligodendrocytes; associated with myelination and neuroprotection.
Key Marker Genes:
MBP (Myelin Basic Protein) 
PLP1 (Proteolipid Protein 1) 
CLDND1 (Intercellular connectivity) 
CRYAB (Heat shock protein) 

Layer 2: Gray Matter / Neurons
Characteristics: Dominated by neurons; associated with nerve signaling and synaptic function.
Key Marker Genes:
SNAP25 (Synaptic vesicle fusion) 
SYT1 (Ca-dependent neurotransmitter release) 
VSNL1 (Neuronal Ca2+ sensor) 

Layer 1: Mixed Region
Characteristics: Represents a transition or mixed region containing both white matter (oligodendrocytes) and gray matter (neurons).

Layer 0: Non-Neuronal / Low Activity
Characteristics: Lacks significant neuronal or myelin-related functions. This region may be rich in other cell types like astrocytes or microglia, or represents less functionally active brain tissue.


🔗 References

Nature Methods, Vol 18 No. 1 (Jan 2021) 
Chitra, U., et al. "Mapping the topography of spatial gene expression with interpretable deep learning." 
Ortiz C et al., Annu Rev Neurosci., 2021 
Mohammadi, E., et al. "Size matters: the impact of nucleus size on results from spatial transcriptomics."
