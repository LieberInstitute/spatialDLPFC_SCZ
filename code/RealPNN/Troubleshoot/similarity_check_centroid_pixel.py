import pandas as pd
import numpy as np
import cv2

# checking if the centroid and pixel counts are different
csv_A1 = '/dcs04/lieber/marmaypag/spatialDLPFC_SCZ_LIBD4100/processed-data/spaceranger/V12F14-053_A1/outs/spatial/tissue_spot_counts.csv'

# Read the dataframe
A1 = pd.read_csv(csv_A1)

# Initialize an empty list to store the column numbers where values are different
different_columns = []

# Iterate over the rows
for index, row in A1.iterrows():
    if row['NDAPI'] != row['CNDAPI']:
        different_columns.append(index)

# Output the column numbers with different values
print(different_columns)
