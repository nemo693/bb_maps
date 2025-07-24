## How to Run the Pipeline

This section provides a step-by-step guide on how to set up the environment and run the entire forest damage detection pipeline.

### 1. Environment Setup

The pipeline requires a specific Conda environment to be activated. The environment is defined in the `fordead_plain.yml` file and can be set up using the provided shell script.

1.  **Navigate to the setup directory:**
    ```bash
    cd setup_environment_fordead
    ```

2.  **Run the installation script:**
    ```bash
    ./install_conda.sh
    ```
    This script will install Miniconda, create the `fordead_plain` environment, and install all the necessary dependencies.

### 2. Pipeline Execution

Once the environment is set up, you can run the main pipeline script.

1.  **Activate the Conda environment:**
    ```bash
    conda activate fordead_plain
    ```

2.  **Navigate to the R scripts directory:**
    ```bash
    cd ../Scripts_fordead_forest_service
    ```

3.  **Run the main orchestrator script:**
    ```bash
    Rscript 00_RUN_Province_Damage_Update.R
    ```
    This will execute the entire pipeline, from data preparation to generating the final damage products.
