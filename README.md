# Glycemic Comparison Index (GCI): A Retrospective Analysis of its Prognostic Value in ICU Patients with AMI and Diabetes

This repository contains the source code for the paper titled "Glycemic Comparison Index (GCI): A Retrospective Analysis of its Prognostic Value in ICU Patients with AMI and Diabetes" (doi:10.1186/s12902-025-01907-2). The code within this repository is designed to replicate the analyses conducted in the study using the MIMIC-IV database.

https://pubmed.ncbi.nlm.nih.gov/40140826/

## Prerequisites

Before you begin, ensure you have met the following requirements:
* You have a basic understanding of SQL and R programming languages.
* You have applied for and obtained access to the MIMIC-IV database. For more information on accessing the MIMIC database, please visit [PhysioNet](https://physionet.org/content/mimiciii/).

## Installation and Setup

1. **MIMIC Database Access and Local Deployment:**
   - Follow the instructions provided by PhysioNet to download and locally deploy the MIMIC-IV database. Ensure the database is correctly set up on your local machine or server where you intend to run the analyses.

2. **Generating Derived Tables:**
   - Clone the [MIMIC-Code repository](https://github.com/MIT-LCP/mimic-code/tree/main/mimic-iv) to your local machine.
   - Navigate to the `concepts` directory within the cloned repository.
   - Execute the SQL scripts to generate derived tables. These tables are essential for the subsequent analyses.

3. **Executing SQL Queries:**
   - In this repository, you will find a series of SQL files numbered sequentially (e.g., `01_first_query.sql`, `02_second_query.sql`, etc.).
   - Run these SQL scripts in the order provided against your local instance of the MIMIC-IV database. These scripts are designed to prepare the data for analysis in R.

4. **Running R Analysis:**
   - Once the SQL queries have been executed and the necessary data is prepared, proceed with the R scripts provided in this repository.
   - The R scripts are responsible for performing the statistical analyses and generating the figures as presented in the paper.

## Usage

Follow the steps outlined above to replicate the analysis. The scripts should be executed in the order they are provided to ensure that all dependencies are correctly resolved.

## Contributing

We welcome contributions to this repository. If you have suggestions to improve the code or analyses, please fork the repository and create a pull request, or open an issue with the tag "enhancement".

## License

Please refer to the LICENSE file for details on the usage and distribution of this code.

## Citation

If you use the code or data provided in this repository in your research, please cite the original paper:

> [Full citation of the paper]

## Acknowledgements

We thank the MIT Laboratory for Computational Physiology and the PhysioNet community for providing access to the MIMIC-IV database.

