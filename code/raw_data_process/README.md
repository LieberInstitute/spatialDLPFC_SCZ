# Software Pre-requisite

## R packages

```{library(here)}
library(tidyverse)
library(cellranger) # Extract excel file
library(glue)
```

# How to start

To run the whole sequencing step,

1.  Set your working directory to the current directory i.e. `cd code/raw_data_process/`

2.  Run `00-main.R` on JHPCE.

# Files:

-   `00-main.R`: master R file to run each task in the preprocessing step
-   `01-create_process_meta.R`: create a meta csv file, `sample_meta_path.csv`,
that contains sample-specific file paths
-   `02-set_fastq_softlink.R`: Set up soft link from data storage to SCZ project
-   `03-prep_for_space_ranger.R`: Prepare and run `spaceranger`

-   Utilit Files:
    - `file_path.R`: a master file contains all path needed for the
    preprocessing step. The validation and creationg of folders are in each
    working script
    - `fun_xxx.R`: helper functions

For more detail explaining each scripts, please see <https://github.com/LieberInstitute/spatialDLPFC_SCZ/issues/15>
