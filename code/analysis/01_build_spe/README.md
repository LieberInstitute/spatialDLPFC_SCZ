Files:
* 01-create_spe_wo_spg.R: Create `spe` object for 63 samples
* 02-create_spe_with_spg.R: Adding SPG information to the `spe` object
* 11_donor_dx.R: Create supplementary table comparing donor balance between diagonsis groups.
* fun_import_dx.R: helper function to import donor information, used in 01-create_spe_wo_spg.R


NOTE:

* Due to the techinical challenges to create excel file on JHPCE, Supplementary tables of donor demographics were created on local machine by running the following line in the Terminal.

  ```
  Rscript 11_donor_dx.R > logs/11_donor_dx.txt 2>&1
  ```


* In the finalized spe object, there are missing values, e.g., `NA` or `nan`, as the mean intensity value (`spg_I` SPG columns in `colData()`) for the four SPG channels. The technical reasons for these missing values is because there are no pixel of the corresponding channel is detected in the spot. Hence, when caclulating the average, the denominator, i.e., the number of pixels, is zero.