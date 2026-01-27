# Sort "ssp_results.csv" by "p" to fix the order scrambled by parallel processing, while preserving formatting.

import pandas as pd

file_name = "data/ssp_results.csv"
df = pd.read_csv(file_name, dtype=str)

df_sorted = df.iloc[pd.to_numeric(df["p"]).sort_values().index]

sorted_file_name = file_name.replace(".csv", "_sorted.csv")
df_sorted.to_csv(sorted_file_name, index=False)

print(f"Sorting complete. Saved to '{sorted_file_name}' with formatting preserved.")
