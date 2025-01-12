# GWASPipe configuration file example

The [config_harmonize_VCF_sumstats.yml](config_harmonize_VCF_sumstats.yml) configuration file is used by GWASPipe to harmonize summary statistics in a specific format. Here's a breakdown of its tasks:

Each task has the following default parameter:
* run: indicates that this step should be performed if set to True, not otherwise.

and the execution order is declared in the ordered dict on top.

**Basic check**: This task checks for basic issues with the input data, such as missing values or incorrect formatting. It provides feedback on whether the input data needs to be corrected before proceeding.

**Infer build**: This step infers the build version of the input data. The `build` refers to the specific version of the genetic database used in the study (e.g., GRCh37 or GRCh38). This information is essential for downstream analysis and interpretation of results.

**Fill data**: In this task, GWASPipe fills missing values in the input data based on predefined rules. For example, it might use a reference distribution to fill missing values for certain variables.

It has the following parameters:

* to_fill: A list of column names to fill with data. In this case, it's set to [ 'MLOG10P', 'P', 'Z' ].
* overwrite: Set to False, which means any existing data in the specified columns will not be overwritten.
* extreme: Set to True, which means that extreme values (e.g., infinity) will be handled accordingly.

**Harmonize**: The harmonization step is where GWASPipe standardizes the input data according to specific formats or protocols. This process ensures that all input data conforms to the required format and can be easily integrated with other datasets.

**Sort alphabetically**: This task sorts the input data in alphabetical order based on a specified column (e.g., study ID). This is useful for maintaining consistency when combining data from multiple studies or experiments.

**Write TSV**: Finally, GWASPipe writes the harmonized and sorted data to a tab-separated values (TSV) file. The `fmt` parameter allows users to specify the desired output format, and other parameters control aspects like compression, float formatting, and more.
