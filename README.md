# Gene Differential Expression Analysis Based on Spatial Transcriptome Data

**Group:** 17th Group
**Date:** April 8th, 2025

---

## 📖 Project Overview
This project investigates the **spatial-level differences in gene expression** within the brain. By utilizing **Spatial Transcriptomics** techniques and reproducing the **GASTON** deep learning method, we successfully mapped the topography of tissue slices and identified distinct functional layers through differential expression analysis.

## 📂 Dataset
* **Source:** Gene Expression Omnibus (GEO)
* **Sample ID:** `GSM5591748` (HE_ON_2339).

## 🛠 Methodology
1.  **Data Preprocessing:** Processing spatial coordinates and gene expression vectors.
2.  **GASTON Reproduction:** We utilized the GASTON method (Generative Adversarial Spatial Transcriptomics Networks) to learn the "topography" of the tissue slice, identifying continuous spatial domains and isodepths.
3.  **Clustering:** Post-GASTON processing involved clustering to unveil distinct brain layers.
4.  **Differential Expression Analysis:** We performed analysis to identify marker genes and determine biological functions for each layer.

---

## 📊 Key Results: Layer Characterization

Our analysis identified four distinct layers with specific gene expression signatures:

### **Layer 3: White Matter / Oligodendrocytes**
* **Biological Function:** Dominated by oligodendrocytes; associated with myelination and neuroprotection.
* **Key Marker Genes:**
  * `MBP`: The important protein of myelination.
  * `PLP1`: The main protein of myelination.
  * `CLDND1`: Intercellular connectivity & barrier function.
  * `CRYAB`: Heat shock protein.

### **Layer 2: Gray Matter / Neurons**
* **Biological Function:** A functional area dominated by neurons; associated with nerve signaling and synaptic function.
* **Key Marker Genes:**
  * `SNAP25`: Key protein for synaptic vesicle fusion.
  * `SYT1`: Ca-dependent neurotransmitter release.
  * `VSNL1`: Neuronal Ca2+ sensor protein family.

### **Layer 1: Mixed Region**
* **Biological Function:** Represents a mixed region containing both white matter (oligodendrocytes) and gray matter (neurons).

### **Layer 0: Non-Neuronal / Low Activity**
* **Biological Function:** Lacks significant neuronal or myelin-related functions.
* **Hypothesis:** May be rich in other cell types (e.g., astrocytes or microglia) or represent less functionally active brain tissue.

---

## 🔗 References
1.  *Nature Methods*, "Method of the Year 2020: Spatially resolved transcriptomics", Vol 18 No. 1 (Jan 2021).
2.  Chitra, U., et al. "Mapping the topography of spatial gene expression with interpretable deep learning".
3.  Ortiz C et al., *Annu Rev Neurosci.*, 2021.
4.  Mohammadi, E., et al. "Size matters: the impact of nucleus size on results from spatial transcriptomics".
