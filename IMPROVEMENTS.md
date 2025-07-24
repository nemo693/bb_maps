## Identified Inconsistencies/Areas for Improvement

During the analysis, several inconsistencies and areas for improvement were noted:

-   **Python Script Locations:** The R script `02_Execute_Core_FORDEAD_Processing.R` refers to Python scripts (`fordead_1.py`, `fordead_2.py`, `fordead_3.py`) located in an external directory (`/mnt/CEPH_PROJECTS/WALDSCHAEDEN/working_folder/fordead_git/forest_outputs`). However, the core logic appears to reside in `fordead_plain/steps/stepX_*.py` within the provided directory structure. Without access to the external directory, it is unclear if the called scripts are wrappers around the `fordead_plain` package or contain separate logic. This should be clarified to ensure the documentation accurately reflects the executed code.

-   **Hardcoded Paths:** This is a pervasive issue across multiple R scripts (`04_Refine_And_Mask_Damage_Products.R`, `05_Style_And_Project_Damage_Maps.R`, `06_Integrate_And_Refine_Damage_Products.R`). Many input and output file paths are hardcoded, making the pipeline inflexible and difficult to manage if the directory structure or data sources change. Examples include:

    -   `province.shp`, `lafis.shp`, `base` raster in `04_Refine_And_Mask_Damage_Products.R`.
    -   `old_monthly` and `old_yearly` paths in `06_Integrate_And_Refine_Damage_Products.R` pointing to a specific date.
    -   QAI input file in `06_Integrate_And_Refine_Damage_Products.R` being hardcoded to a single tile and date.

-   **Redundant Logic:**

    -   The calculation of `last_detectable_date` is duplicated in `04_Refine_And_Mask_Damage_Products.R` and `05_Style_And_Project_Damage_Maps.R`. This logic should be centralized.
    -   QAI processing logic (converting bit values, filtering) is present in both `02_Execute_Core_FORDEAD_Processing.R` (R) and `06_Integrate_And_Refine_Damage_Products.R` (R). This redundancy could lead to inconsistencies if not carefully managed.

-   **Parameter Management:** The comment in `06_Integrate_And_Refine_Damage_Products.R` (`# todo: write parameters separately as a file (e.g. outfold keeps reverting to an older setting)`) indicates a known issue with parameter persistence. While Python scripts use `param.py`, a more robust and centralized system for managing parameters across all R scripts would improve maintainability and prevent unexpected behavior.

-   **Fragile Date/Naming Logic:**

    -   The "fix november" and "fix april/winter" adjustments in `04_Refine_And_Mask_Damage_Products.R` and `05_Style_And_Project_Damage_Maps.R` suggest specific workarounds that might indicate underlying issues with date handling or data consistency.
    -   Hardcoded years in output filenames (e.g., `_2025_25832`) in `05_Style_And_Project_Damage_Maps.R` and `06_Integrate_And_Refine_Damage_Products.R` require manual updates for each new year.

-   **Implicit Dependencies:** Several scripts rely on variables defined in previously sourced scripts (e.g., `index.list` in `03_Mosaic_FORDEAD_Outputs.R` from `02_Execute_Core_FORDEAD_Processing.R`), which can make it harder to understand individual script functionality in isolation.

-   **Hardcoded Visuals:** Color palettes and labels are directly embedded in `05_Style_And_Project_Damage_Maps.R`, making them less dynamic and harder to update without modifying the script itself.
