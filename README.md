# Inslulin-mapping

To replicate the computational analysis described in the paper, including phosphoproteomics, metabolomics, and kinetic modeling of insulin signaling pathways, the following R code outlines a comprehensive approach. We will go step by step, covering how to analyze both phosphoproteomics and metabolomics data, reconstruct signaling pathways, and perform dynamic modeling. This code assumes that you have the datasets of phosphoproteomics and metabolomics (e.g., time-course data for insulin-stimulated conditions).

### **Step 1: Load Necessary Libraries**
First, we need to load essential R libraries for data analysis, modeling, and visualization.

```r
# Load necessary libraries
install.packages(c("tidyverse", "limma", "gplots", "igraph", "deSolve", "phosphoR", "KEGGREST"))
library(tidyverse)
library(limma)
library(gplots)
library(igraph)
library(deSolve) # For kinetic modeling
library(phosphoR) # Phosphoproteomics analysis
library(KEGGREST) # For pathway reconstruction
```

### **Step 2: Load Phosphoproteomics and Metabolomics Data**
Load the time-course data for both phosphoproteomics and metabolomics (assume CSV format). We will also preprocess the data (log transformation and normalization).

```r
# Load data (replace with actual file paths)
phospho_data <- read.csv("phosphoproteomics_data.csv")
metabolomics_data <- read.csv("metabolomics_data.csv")

# Log2 transformation for phosphoproteomics
phospho_data_log <- phospho_data %>%
  mutate(across(where(is.numeric), log2))

# Normalization (using limma package)
phospho_data_norm <- normalizeBetweenArrays(as.matrix(phospho_data_log), method = "quantile")

# For metabolomics, normalization
metabolomics_data_log <- metabolomics_data %>%
  mutate(across(where(is.numeric), log2))

# Example output of normalized data
head(phospho_data_norm)
head(metabolomics_data_log)
```

### **Step 3: Differential Phosphorylation Analysis**
Identify significantly phosphorylated proteins at different time points after insulin stimulation.

```r
# Create design matrix for time-course experiment (adjust based on your study design)
time_points <- factor(c("T0", "T10", "T30", "T60")) # T0 as control, others as insulin treatment
design_matrix <- model.matrix(~ time_points)

# Fit the linear model
fit <- lmFit(phospho_data_norm, design_matrix)

# Calculate differential expression (phosphorylation) for each time point compared to T0
fit <- eBayes(fit)
results <- topTable(fit, coef = 2:4, adjust.method = "BH", number = Inf)

# Filter significant proteins (adjusted p-value < 0.05)
significant_proteins <- results %>%
  filter(adj.P.Val < 0.05)

# Output significant results
head(significant_proteins)
```

### **Step 4: Metabolomics Analysis**
Identify significant changes in metabolite levels.

```r
# Create design matrix for metabolomics time-course data
fit_met <- lmFit(metabolomics_data_log, design_matrix)
fit_met <- eBayes(fit_met)
met_results <- topTable(fit_met, coef = 2:4, adjust.method = "BH", number = Inf)

# Filter significant metabolites (adjusted p-value < 0.05)
significant_metabolites <- met_results %>%
  filter(adj.P.Val < 0.05)

# Output significant metabolites
head(significant_metabolites)
```

### **Step 5: Pathway Enrichment Analysis**
Identify pathways enriched by the differentially phosphorylated proteins and metabolites using the KEGG database.

```r
# Use KEGGREST to retrieve KEGG pathway data
# Map significant proteins to KEGG pathway identifiers
kegg_protein_ids <- keggConv("hsa", significant_proteins$Gene)
kegg_metabolite_ids <- keggConv("hsa", significant_metabolites$Metabolite)

# Pathway enrichment analysis for proteins
kegg_pathways_protein <- keggList("pathway", "hsa")
pathway_hits_protein <- sapply(kegg_protein_ids, function(id) keggLink("pathway", id))

# Pathway enrichment analysis for metabolites
kegg_pathways_metabolite <- keggList("pathway", "hsa")
pathway_hits_metabolite <- sapply(kegg_metabolite_ids, function(id) keggLink("pathway", id))

# Output top pathways for both proteins and metabolites
head(pathway_hits_protein)
head(pathway_hits_metabolite)
```

### **Step 6: Reconstruct Signal Flow**
Reconstruct the insulin signal flow network based on the KEGG pathways for significant proteins and metabolites.

```r
# Build the network of proteins and metabolites involved in insulin signaling
# Create an igraph object for the network
nodes <- data.frame(name = c(significant_proteins$Gene, significant_metabolites$Metabolite))
edges <- data.frame(from = c(significant_proteins$Gene), to = c(significant_metabolites$Metabolite))

# Create graph object
g <- graph_from_data_frame(edges, vertices = nodes, directed = TRUE)

# Plot the network
plot(g, vertex.label = V(g)$name, vertex.size = 5, edge.arrow.size = 0.5, layout = layout.fruchterman.reingold)
```

### **Step 7: Dynamic Modeling of Glycolytic Pathway**
Develop a kinetic model to simulate the dynamic signal flow in glycolysis under insulin stimulation. This example focuses on the conversion of fructose-6-phosphate (F6P) to fructose-1,6-bisphosphate (F1,6BP).

```r
# Define differential equations for glycolytic pathway
glycolysis_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Equations for F6P and F1,6BP
    dF6P <- -k1 * F6P + k2 * F1,6BP
    dF1_6BP <- k1 * F6P - k2 * F1,6BP
    
    # Return the rate of change
    return(list(c(dF6P, dF1_6BP)))
  })
}

# Set initial conditions and parameters
initial_state <- c(F6P = 1, F1_6BP = 0)
parameters <- c(k1 = 0.5, k2 = 0.1)

# Simulate the system over time
time <- seq(0, 60, by = 1)
output <- ode(y = initial_state, times = time, func = glycolysis_model, parms = parameters)

# Plot the results
output_df <- as.data.frame(output)
ggplot(output_df, aes(x = time)) +
  geom_line(aes(y = F6P, color = "F6P")) +
  geom_line(aes(y = F1_6BP, color = "F1,6BP")) +
  labs(title = "Dynamic Simulation of Glycolysis", y = "Concentration", x = "Time (minutes)") +
  theme_minimal()
```

### **Step 8: Feedback Mechanisms and Allosteric Regulation**
Add feedback loops and allosteric regulation based on metabolites like PFKL and others.

```r
# Modify the kinetic model to include feedback from allosteric effectors
glycolysis_feedback_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Allosteric feedback from downstream metabolites (e.g., PEP, citrate)
    feedback_PEP <- 1 / (1 + PEP / K_PEP)
    feedback_citrate <- 1 / (1 + citrate / K_citrate)
    
    # Adjust reaction rates with feedback
    dF6P <- -k1 * F6P * feedback_PEP * feedback_citrate + k2 * F1,6BP
    dF1_6BP <- k1 * F6P * feedback_PEP * feedback_citrate - k2 * F1,6BP
    
    return(list(c(dF6P, dF1_6BP)))
  })
}

# Simulate feedback system
output_feedback <- ode(y = initial_state, times = time, func = glycolysis_feedback_model, parms = parameters)

# Plot results
output_feedback_df <- as.data.frame(output_feedback)
ggplot(output_feedback_df, aes(x = time)) +
  geom_line(aes(y = F6P, color = "F6P")) +
  geom_line(aes(y = F1_6BP, color = "F1,6BP")) +
  labs(title = "Glycolysis with Feedback Regulation", y = "Concentration", x = "Time (minutes)") +
  theme_minimal()
```

### **Conclusion**
This comprehensive R code allows for the analysis of phosphoproteomics and metabolomics data, pathway enrichment, signal flow reconstruction, and dynamic modeling of glycolysis under insulin signaling. You can modify and extend the code based on the specifics of your datasets and research questions.
