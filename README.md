# HIV-1-Molecular-Explorer
Molecular epidemiology/Sequence retrieval/Multiple alignment/Phylogenetic analysis

## Overview
The HIV-1 Molecular Explorer is an interactive application designed to search, collect, analyze, and visualize HIV-1 genetic data directly from the National Center for Biotechnology Information (NCBI) database. Built for researchers and bioinformatics enthusiasts, this MVP simplifies the process of data acquisition, sequence alignment, and evolutionary visualization.

## Key Features
- **Targeted Data Retrieval:** Fetches HIV-1-related data specifically from the past 10 years using the NCBI API.
- **Customizable Search Parameters:** Allows users to filter data by sequence length, target genes, country of origin, and year.
- **Comprehensive Data Export:** Automatically parses and saves collected datasets in CSV format. The output includes:
  - Year of Isolation
  - Title
  - Genbank ID
  - Country
  - Viral DNA Sequence
- **Advanced Sequence Alignment:** Features a built-in sequence alignment tool.
- **Phylogenetic Visualization:** Easily export sequence alignment results as a phylogenetic tree.

## Setup & Configuration
**Important NCBI API Requirement:** 
To comply with NCBI's usage policies and prevent your IP/API calls from being blocked, Entrez requires a registered email string.
*   **Action Required:** Ensure you enter a valid email into the designated sidebar input within the app before running a search.

## Verification & Testing
*   **Data Integrity Check:** Built-in test script to automatically verify that the fetched data is intact, formatted correctly, and that the core functions (API fetching, CSV saving, alignment) execute correctly without silent errors.

## Usage Guide
1. Launch the application and navigate to the sidebar.
2. Enter your email address.
3. Input your desired search parameters.
4. Execute the search to collect data and review the generated dataset.
5. Save the raw dataset as a CSV or proceed to the alignment tab.
6. Run the sequence alignment and export your resulting phylogenetic tree.
