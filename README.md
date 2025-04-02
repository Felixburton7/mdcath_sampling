flowchart TD
    A["Molecular Dynamics Data<br>CSVs at 5 Temps + Avg"] 
    B["CATH Classifications v4.4.0"]
    
    A & B --> D["Extract CATH Hierarchy<br>Class, Architecture, Topology, Homology"]
    A --> E["Aggregate Domain-Level Features"]
    E --> E1["Secondary Structure %<br>Helix, Sheet, Coil"]
    E --> E2["Backbone Flexibility<br>Average RMSF"]
    E --> E3["Core/Exterior Ratio"]
    E --> E4["Solvent Accessibility"]
    E --> F["Temperature Stability Profiles<br>Stable: RMSF450/RMSF320 < 1.0<br>Moderate: 1.0-2.0<br>Unstable: > 2.0"]
    
    E1 & E2 & E3 & E4 & F --> G["Hierarchical Stratification by CATH"]
    D --> G
    G --> H["Targeted Sampling within Topology Groups"]
    H --> H1{"Topology Size ≥ 10 domains?"}
    
    H1 -->|Yes| H2["K-means Clustering on Feature Vectors"]
    H2 --> H3["Select One Representative Domain per Cluster"]
    
    H1 -->|No| H4["Proportional Sampling Based on Feature Uniqueness"]
    
    H3 & H4 --> I["Homology Control"]
    I --> I1["Build Domain Network<br>PDB Origin / Homologous Superfamily"]
    I1 --> I2["Sample Connected Components<br>Keep Related Domains Together"]
    
    I2 --> J["Statistical Validation"]
    J --> J1["KS Tests for Continuous Distributions<br>Domain Size, RMSF Values"]
    J --> J2["χ² Tests for Categorical Distributions<br>CATH Classes, Secondary Structure"]
    J1 & J2 --> K["Calculate Representation Index RI<br>Geometric Mean of Similarity Metrics"]
    K --> L{"RI ≥ 0.9?"}
    L -->|No| M["Refine Sampling"]
    M --> G
    
    L -->|Yes| N["Holdout Set<br>10%, 540 Domains"]
    L -->|Yes| O["Training Set<br>90%, 4,858 Domains"]
    N & O --> P["Output Domain Lists to Text Files"]
    
    classDef input fill:#D4F1F9,stroke:#05386B,stroke-width:2px,rx:5px;
    classDef process fill:#E0F4FF,stroke:#05386B,stroke-width:2px,rx:5px;
    classDef feature fill:#CCE8FF,stroke:#05386B,stroke-width:1px,rx:5px;
    classDef decision fill:#FFEECC,stroke:#D9730D,stroke-width:2px,rx:10px;
    classDef validation fill:#F0E6FF,stroke:#5E2CA5,stroke-width:2px,rx:5px;
    classDef output fill:#E3FCEF,stroke:#0A7B83,stroke-width:2px,rx:5px;
    
    class A,B input;
    class D,E,G,H,H2,H3,H4,I,I1,I2,M process;
    class E1,E2,E3,E4,F feature;
    class H1,L decision;
    class J,J1,J2,K validation;
    class N,O,P output;
