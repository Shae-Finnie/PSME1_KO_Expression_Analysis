#!/usr/bin/env python3
#
# 00_csv_prep.py
# Extracts PSME1 and scram data from original unnormalized expression
# csv
#
# Usage: python scripts/00_csv_prep.py
# run from PSME1_KO_Expression_Analysis root 

import pandas as pd
# Load specific cols by index
PSME1_vs_Scram = pd.read_csv(
    "data/raw/PSME1_vs_scram_raw.csv", usecols=[0,4,5,6,7,8,9]
    )
# Save new csv
PSME1_vs_Scram.to_csv("data/raw/PSME1_vs_Scram.csv",index=False)
