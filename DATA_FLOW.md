## Data Flow

The data flows through the application in a sequential manner:

1.  **Raw Satellite Data (Sentinel-2 BOA & QAI):** Stored in `in_path` (e.g., `/mnt/CEPH_PROJECTS/sao/SENTINEL-2/SentinelVegetationProducts/FORCE/level2/`).
2.  **`01_Import_S2_data.R`:** Organizes raw data into grid-specific `_boa` and `_update` directories.
3.  **`02_Execute_Core_FORDEAD_Processing.R`:**
    -   Reads organized BOA and QAI data.
    -   Generates R-based QAI masks.
    -   Passes parameters to Python scripts via `param.py`.
    -   Python wrapper scripts (`fordead_1.py`, `fordead_2.py`, `fordead_3.py`) process data, generating:
        -   Vegetation indices and masks (`fordead_output_*/VegetationIndex`, `fordead_output_*/Mask`).
        -   Model coefficients (`fordead_output_*/DataModel/coeff_model.tif`).
        -   Anomaly, dieback, and stress data (`fordead_output_*/DataAnomalies`, `fordead_output_*/DataDieback`, `fordead_output_*/DataStress`).
4.  **`03_Mosaic_FORDEAD_Outputs.R`:** Reads tiled outputs from the Python steps and mosaics them into full-area raster and vector files (e.g., `output_NDWI_merged.shp`, `output_confidence_NDVI_merged.tif`).
5.  **`04_Refine_And_Mask_Damage_Products.R`:** Reads merged outputs, applies further masking (NDWI, province, LAFIS), and generates refined monthly and annual products (e.g., `NDVI_merged_masked_with_NDWI_province_lafis.tif`).
6.  **`05_Style_And_Project_Damage_Maps.R`:** Reads the refined monthly and annual products, applies styling, and generates initial final GeoTIFFs and Shapefiles (e.g., `changes_yearly_damages_july_2025_25832.tif`).
7.  **`06_Integrate_And_Refine_Damage_Products.R`:** Reads the initial final products, integrates historical data, applies additional filters (stress, shadow, single pixel), and produces the definitive final GeoTIFFs and Shapefiles (e.g., `changes_monthly_damages_july_2025_25832_final.tif`).
8.  **`00_RUN_Province_Damage_Update.R`:** Uses the final GeoTIFFs for `gdaladdo` and `scp` commands to generate overviews and upload data.

## Data Flow Diagram

``` mermaid
graph TD
    subgraph "R Scripts"
        A["00_RUN_Province_Damage_Update.R (Orchestrator)"] --> B["01_Import_S2_data.R"];
        B --> C["02_Execute_Core_FORDEAD_Processing.R"];
        C --> D{Python Scripts};
        D --> E["03_Mosaic_FORDEAD_Outputs.R"];
        E --> F["04_Refine_And_Mask_Damage_Products.R"];
        F --> G["05_Style_And_Project_Damage_Maps.R"];
        G --> H["06_Integrate_And_Refine_Damage_Products.R"];
    end

    subgraph "Python Core (fordead_plain)"
        D -- calls --> Py1["fordead_1.py (wrapper)"];
        Py1 -- executes --> Step1["step1_compute_masked_vegetationindex.py"];
        D -- calls --> Py2["fordead_2.py (wrapper)"];
        Py2 -- executes --> Step2["step2_train_model.py"];
        D -- calls --> Py3["fordead_3.py (wrapper)"];
        Py3 -- executes --> Step3["step3_dieback_detection.py"];
    end

    subgraph "Inputs"
        I[Raw Sentinel-2 Data] --> B;
        J[Previous Run Outputs] --> A;
    end

    subgraph "Outputs"
        H --> K["Final Damage Products (GeoTIFF, Shapefile)"];
    end

    style D fill:#f9f,stroke:#333,stroke-width:2px
```
