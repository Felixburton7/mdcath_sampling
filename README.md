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
