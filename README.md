# **Gene Expression Analysis of Influenza A’s Impact on Plasmacytoid Dendritic Cells (pDCs)**

This repository contains the code and analysis pipeline for processing gene expression data from the **Gene Expression Omnibus (GEO)** dataset **GSE68849**. The dataset includes information on the effects of **Influenza A** infection on **plasmacytoid dendritic cells (pDCs)**. The goal of the analysis is to identify **differentially expressed genes (DEGs)**—genes whose expression levels change significantly in response to the virus.

---

## **Project Overview**

The analysis focuses on comparing **control samples** (uninfected pDCs) with **treated samples** (pDCs infected with Influenza A). This allows us to identify genes whose expression is significantly altered by the infection, which may provide insights into the molecular mechanisms underlying the immune response to the virus.

---

## **Key Steps in the Analysis**

### **Step 1: Downloading, Extracting, and Inspecting Data**

- The raw gene expression data is provided as a compressed `.gz` file. The first step is to decompress the file and inspect the structure to extract the relevant expression matrix.
  
- The gene expression matrix starts at line 61 of the file, so we clean the data by extracting the relevant portion and removing metadata.

### **Step 2: Data Cleaning (Bash)**

- Cleaning up unnecessary characters, such as double quotes in the data, ensures that the dataset is ready for analysis.

### **Step 3: Setting Up Python for Differential Expression Analysis**

- We use **Python** and the **pandas** and **scipy** libraries to load the cleaned data, prepare it for statistical analysis, and define the necessary treatment groups (**control vs. treated**).

- The metadata is used to assign samples to the control and treated groups.

### **Step 4: Assigning Control and Treated Samples**

- We extract the metadata from the dataset to map the samples to their respective **control** and **treated** categories.

### **Step 5: Data Validation and Cleaning**

- After assigning control and treated groups, we validate the data by removing extra quotes from sample names and filtering out only the control and treated samples for the analysis.

### **Step 6: Statistical Analysis (Differential Expression)**

- We perform a **t-test** to identify differentially expressed genes (DEGs) by comparing the expression levels between the control and treated groups.

- The results are corrected for multiple testing using the **Benjamini-Hochberg** procedure to control the **false discovery rate (FDR)**.

### **Step 7: Visualization**

- A **volcano plot** is generated to visualize the fold change and significance of DEGs.

- A **heatmap** is created to visualize the expression patterns of the top 20 DEGs across samples.

### **Step 8: Pathway Enrichment Analysis**

- The **GProfiler** package is used to identify biological pathways enriched among the DEGs, providing a biological context for the gene expression changes observed.

---

## **Dependencies**

To run the analysis, you need the following Python packages:

- **pandas**
- **scipy**
- **seaborn**
- **matplotlib**
- **gprofiler**

---

## **Data Input and Output**

- **Input Data**: A gene expression matrix from the **GSE68849** dataset (typically in `.txt` or `.gz` format).

- **Output Data**:
  - A CSV file (**differentially_expressed_genes.csv**) containing the differentially expressed genes.
  - Visualization files such as **volcano plots** and **heatmaps** showing the results of the analysis.
  - Pathway enrichment results from **GProfiler**.

---

## **Conclusion**

This analysis pipeline processes gene expression data to identify how **Influenza A** infection impacts gene expression in **plasmacytoid dendritic cells**. The results provide valuable insights into the immune response to the virus, with potential applications in understanding the molecular mechanisms of influenza pathogenesis.
