# mdCATH Dataset Stratified Sampling

This repository implements a hierarchical stratified sampling protocol to partition the mdCATH dataset into training (90%) and holdout (10%) sets. The goal is to preserve the diversity of structural and dynamic properties across protein domains while ensuring that homologous domains remain together.

---

## Overview

The mdCATH dataset is preprocessed by:
- **Extracting CATH Classifications**: Domains are classified by class, architecture, topology, and homology using CATH v4.4.0.
- **Feature Aggregation**: Domain-level features are computed from residue-level molecular dynamics simulations at five different temperatures (320K–450K) plus an average. Key features include:
  - **Secondary Structure Composition**: Percentages of helix, sheet, and coil.
  - **Backbone Flexibility**: Average RMSF values.
  - **Core/Exterior Ratio**: Proportion of residues in the core.
  - **Solvent Accessibility**: Average relative accessibility.
- **Temperature Stability Profiles**: Domains are categorized as:
  - **Stable**: RMSF₄₅₀/RMSF₃₂₀ < 1.0
  - **Moderately Stable**: Ratio between 1.0 and 2.0
  - **Unstable**: Ratio > 2.0

---

## Stratified Sampling Protocol

### Hierarchical Stratification
- **Partitioning**: Domains are first grouped based on the CATH classification hierarchy.
- **Targeted Sampling within Topology Groups**:
  - For topologies with ≥10 domains:  
    - **K-means clustering** is applied to the feature vectors.
    - One representative domain per cluster is selected to maximize diversity.
  - For smaller topologies:
    - Domains are sampled with a probability proportional to their uniqueness in feature space.

### Homology Control
- **Domain Network Construction**:  
  - A network is built where nodes (domains) are connected if they share the same PDB origin or belong to the same homologous superfamily.
- **Component Sampling**:  
  - Instead of sampling individual domains, connected components are sampled to ensure that related domains remain together in either the training or holdout set.

### Validation and Iterative Refinement
- **Statistical Validation**:
  - Kolmogorov-Smirnov tests for continuous distributions (e.g., domain size, RMSF values).
  - χ² tests for categorical distributions (e.g., CATH classifications, secondary structure composition).
- **Representation Index (RI)**:
  - A composite RI (geometric mean of distribution similarity, hierarchy coverage, and stability profile coverage) is computed.
  - The sampling process is iteratively refined until RI ≥ 0.9.
- **Final Outcome**:
  - **Holdout Set**: 540 domains (10% of the dataset).
  - **Training Set**: 4,858 domains.

---

## Code Overview

The repository includes several key Python modules:

### Data Loading & Feature Aggregation
- **`load_temperature_data(file_paths)`**:  
  Loads CSV files containing molecular dynamics simulation data for each temperature.
- **`aggregate_domain_features(temperature_data)`**:  
  Aggregates residue-level data into domain-level features including secondary structure, RMSF statistics, core/exterior ratio, and solvent accessibility. It also computes temperature stability profiles.

### Hierarchical Stratified Sampling
- **`hierarchical_stratified_sampling(domain_features, cath_dict, sample_ratio=0.1)`**:  
  Implements the multi-level stratified sampling by grouping domains based on the CATH hierarchy, applying k-means clustering for larger groups, and sampling proportionally for smaller groups.
- **`network_aware_sampling(candidate_domains, domain_features, cath_dict, sample_ratio=0.1)`**:  
  Applies network-based homology control by constructing a domain network and sampling connected components.

### Main Pipeline
- The main script orchestrates the process:
  1. **Parse CATH Classifications**: Reads the CATH file.
  2. **Load Temperature Data**: Reads simulation data from CSV files.
  3. **Compute Domain-Level Features**: Aggregates features from the loaded data.
  4. **Perform Stratified Sampling**: Applies hierarchical and network-aware sampling.
  5. **Validation and Refinement**: Validates the holdout set using statistical tests and refines the sampling iteratively until a robust Representation Index (RI) is achieved.
  6. **Output**: Saves the holdout and training domain IDs into separate text files.

### Example Execution
```python
import os
import time
from cath_parser import parse_cath_file
from feature_extractor import load_temperature_data, aggregate_domain_features
from sampling import hierarchical_stratified_sampling
from validation import statistical_validation

def main():
    print("Starting mdCATH Dataset Stratified Sampling...")
    
    # Parse CATH classifications
    cath_dict = parse_cath_file('/path/to/cath-domain-list.txt')
    
    # Load temperature data
    file_paths = [
        'final_dataset_temperature_320.csv',
        'final_dataset_temperature_348.csv',
        'final_dataset_temperature_379.csv',
        'final_dataset_temperature_413.csv',
        'final_dataset_temperature_450.csv',
        'final_dataset_temperature_average.csv'
    ]
    file_paths = [os.path.join('/path/to/ML_features', f) for f in file_paths]
    temperature_data = load_temperature_data(file_paths)
    
    # Compute domain-level features
    domain_features = aggregate_domain_features(temperature_data)
    
    # Perform hierarchical stratified sampling
    holdout_domains = hierarchical_stratified_sampling(domain_features, cath_dict, sample_ratio=0.1)
    
    # Validate the holdout set
    all_domains = list(domain_features.keys())
    validation_results = statistical_validation(holdout_domains, all_domains, domain_features, cath_dict)
    ri = validation_results['representation_index']
    
    # Iterative refinement if necessary...
    
    # Save holdout and training domains
    with open('mdcath_holdout_domains.txt', 'w') as f:
        for domain in sorted(holdout_domains):
            f.write(f"{domain}\n")
    
    training_domains = list(set(all_domains) - set(holdout_domains))
    with open('mdcath_training_domains.txt', 'w') as f:
        for domain in sorted(training_domains):
            f.write(f"{domain}\n")
    
    print("Stratified sampling completed successfully!")

if __name__ == "__main__":
    start_time = time.time()
    main()
    elapsed_time = time.time() - start_time
    print(f"Total execution time: {elapsed_time:.2f} seconds")






MEDUSA
graph TD
    %% Input Data
    input[MD Simulations at 5 Temperatures + CATH v4.4.0]
    
    %% Feature Extraction
    input --> features[Domain-Level Feature Extraction]
    features --> feat1[Secondary Structure %]
    features --> feat2[Backbone Flexibility (RMSF)]
    features --> feat3[Core/Exterior Ratio]
    features --> feat4[Solvent Accessibility]
    features --> stability[Temperature Stability Profiles]
    
    %% Core Sampling Process
    features --> strat[Hierarchical Stratification by CATH]
    strat --> topology{Topology Size}
    
    %% Large Topology Groups
    topology -->|≥10 domains| kmeans[K-means Clustering]
    kmeans --> rep[Select Representative Domains]
    
    %% Small Topology Groups
    topology -->|<10 domains| prop[Proportional Sampling by Uniqueness]
    
    %% Homology Control
    rep --> homology[Homology Control]
    prop --> homology
    homology --> network[Build Domain Network]
    network --> components[Sample Connected Components]
    
    %% Validation Loop
    components --> validate[Statistical Validation]
    validate --> ri[Representation Index (RI)]
    ri --> threshold{RI ≥ 0.9?}
    threshold -->|No| refine[Refine Sampling]
    refine --> strat
    
    %% Final Output
    threshold -->|Yes| output[10% Holdout Set (540 domains)]
    
    %% Add Visual Clarity
    classDef process fill:#D4F1F9,stroke:#05386B,stroke-width:2px
    classDef decision fill:#FFEECC,stroke:#D9730D,stroke-width:2px
    classDef output fill:#E3FCEF,stroke:#0A7B83,stroke-width:2px
    
    class input,features,feat1,feat2,feat3,feat4,stability,strat,kmeans,rep,prop,homology,network,components,validate,ri,refine process
    class topology,threshold decision
    class output output

graph TD
    subgraph Inputs
        direction LR
        A[Molecular Dynamics Data (CSVs @ 5 Temps + Avg)]
        B[CATH Classifications (v4.4.0)]
        C[Raw mdCATH Domains (Implicit)]
    end

    subgraph Preprocessing & Feature Engineering
        direction TB
        D[Extract CATH Hierarchy (C, A, T, H)]
        E[Aggregate Domain-Level Features]
        E --> E1(Sec. Structure %)
        E --> E2(Avg. RMSF)
        E --> E3(Core/Exterior Ratio)
        E --> E4(Avg. Rel. Solvent Acc.)
        F[Calculate Temperature Stability Profile (Stable/Moderate/Unstable)]
    end

    subgraph Core Sampling Protocol
        direction TB
        G[Hierarchical Stratification by CATH]
        H[Targeted Sampling within Topology Groups]
        H --> H1{>= 10 Domains?}
        H1 -- Yes --> H2[K-means Clustering on Features]
        H2 --> H3[Select 1 Representative/Cluster]
        H1 -- No --> H4[Proportional Sampling (Uniqueness)]
        I[Homology Control]
        I --> I1[Build Domain Network (PDB Origin / Homology)]
        I1 --> I2[Sample Connected Components]
    end

    subgraph Validation & Refinement Loop
        direction TB
        J[Statistical Validation]
        J --> J1(Kolmogorov-Smirnov Tests - Continuous)
        J --> J2(Chi-squared Tests - Categorical)
        K[Calculate Representation Index (RI)]
        K --> K1(Composite: Dist. Sim., Hier. Cov., Stab. Cov.)
        L{RI >= 0.9?}
        L -- No --> CoreSampling[Refine Sampling - Back to Core Protocol]
        L -- Yes --> Outputs[Proceed to Final Output]
    end

    subgraph Final Outputs
        direction LR
        M[Holdout Set (10%, 540 Domains)]
        N[Training Set (90%, 4858 Domains)]
        O[Output Files (`holdout_domains.txt`, `training_domains.txt`)]
    end

    subgraph Implementation Details
        direction TB
        P[Python Code Modules]
        P --> P1(`load_temperature_data`)
        P --> P2(`aggregate_domain_features`)
        P --> P3(`hierarchical_stratified_sampling`)
        P --> P4(`network_aware_sampling`)
        P --> P5(`statistical_validation`)
        Q[Main Pipeline Script (Orchestration)]
    end

    %% Central Process Hub
    Hub((mdCATH Stratified Sampling))

    %% Connections
    A --> Hub
    B --> Hub
    C --> Hub

    Hub --> D
    Hub --> E
    D --> G
    E --> H
    E --> F
    F --> K1  % Stability profile feeds into RI calculation

    Hub --> G
    G --> H
    H --> I2 % Output of targeted sampling feeds into component sampling pool
    Hub --> I  % Homology control is a core part

    I2 --> J % Sampled set goes to validation
    K --> L
    J --> K

    L -- Yes --> M
    L -- Yes --> N
    M --> O
    N --> O

    %% Implementation supports the process
    Implementation_Details --> Hub

    %% Linking Preprocessing Outputs to Usage
    D --> G  % CATH hierarchy used for stratification
    E --> H  % Features used for clustering/proportional sampling
    E --> J  % Features used for statistical validation distributions
    D --> J  % CATH classifications used for categorical validation

    %% Link Core Sampling output to Validation
    I2 --> J % The result of component sampling is validated
