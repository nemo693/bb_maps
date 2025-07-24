# Forest Damage Detection Pipeline

This project contains a multi-stage pipeline for detecting forest damage and mortality using Sentinel-2 satellite imagery. The core processing logic is implemented in a Python package (`fordead_plain`), while the overall orchestration, data preparation, and post-processing are handled by a series of R scripts.

This documentation is intended to guide a new user through the pipeline's architecture, setup, and execution.

## Documentation

For a complete understanding of the pipeline, please refer to the following documents:

- **[OVERVIEW.md](./OVERVIEW.md):** A high-level summary of the application's purpose and design.
- **[HOW_TO_RUN.md](./HOW_TO_RUN.md):** Step-by-step instructions for setting up the environment and executing the pipeline.
- **[DATA_FLOW.md](./DATA_FLOW.md):** A detailed description and diagram of the data flow between the different scripts and processes.
- **[R_SCRIPTS.md](./R_SCRIPTS.md):** Detailed documentation for each of the R scripts that control the pipeline.
- **[PYTHON_PACKAGE.md](./PYTHON_PACKAGE.md):** Detailed documentation for the core Python package and its wrapper scripts.
- **[IMPROVEMENTS.md](./IMPROVEMENTS.md):** A list of identified inconsistencies and areas for potential improvement in the codebase.
