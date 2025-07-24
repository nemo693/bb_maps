# This README describes two modes of operation (A and B) for the forest damage detection process.

# Mode A: Incremental Update
# This mode runs within a specific conda environment where the core code has been modified.
# The model is trained on historical data from year y(0) up to y(t-1).
# It is then used to detect changes occurring specifically during year y(t).
# This mode is designed to work only with mask versions that have already been processed.
A: runs in a conda environment where the codes have been modified (see code in the folder). The model is trained for years y(0) to y(t-1) and it is used to detect changes occurring during y(t). Only works with mask versions which have already been processed.

# Mode B: Full Retraining and Detection
# This mode also runs within a modified conda environment.
# The model is trained on a broader historical dataset, from year y(0) up to year y(n).
# It is then used to detect changes occurring in subsequent years, from y(n) up to y(t).
# This mode can also be used to generate new mask versions, implying it's suitable for initial setup or major retraining.
B: runs in a conda environment where the codes have been modified (see code in the folder). The model is trained for years y(0) to y(n) and it is used to detect changes occurring during the following years y(n) - y(t). Can also be used to generate mask versions. 


