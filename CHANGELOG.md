# MetaflowX: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

# Changelog for MetaflowX v0.1.0

## v0.1.0 - [2024-11-05]
### Initial Release
We are excited to announce the first official release of MetaflowX! MetaflowX is an innovative workflow management tool that merges the power of Nextflow with the flexibility of Python, offering efficient and scalable solutions for data pipelines. Below are the key features and highlights of version v1.0.0:

### Features
- **Workflow Definition and Management**: Enables users to define complex data workflows using straightforward Python code.
- **Nextflow Integration**: Deep integration with Nextflow for automated task parallelization and resource management.
- **Scalability**: Supports distributed execution, allowing seamless scaling in both local and cloud environments.
- **Intuitive Debugging Tools**: Provides detailed logs and error messages for quick troubleshooting.
- **Automated Data Tracking**: Built-in data versioning to track inputs and outputs for each task.
- **Flexible Parameter Management**: Custom parameter injection for simplified configuration across various runtime scenarios.
- **Comprehensive Monitoring Dashboard**: Real-time monitoring of workflow status, resource utilization, and performance metrics.

### Technical Highlights
- Compatible with Python 3.8 and newer versions to ensure optimal performance in modern Python environments.
- Built-in support for cloud integrations such as AWS S3 and GCP Storage for seamless data access and persistence.
- Modular architecture designed for future expansion and easy integration of new components and plugins.

### Known Issues
- High concurrency scenarios may require further optimization for resource scheduling.
- Log output might experience slight delays under certain edge cases.

---

For more detailed documentation and usage guides, please visit the 🚀 [MetaflowX User Manual](README.md).

# Changelog for MetaflowX v1.0.0

## v1.0.0 - [2025-08-04]
We are excited to announce the first official release of MetaflowX v1.0.0! This update brings major improvements and new features to the MetaflowX ecosystem.

Key updates include:

- Added a stub test and demo pipeline for simplified validation and quick-start testing
- Updated the PBO and SBO modules with improved performance and compatibility
- Enhanced support for multi-threaded CD-HIT, enabling faster and more scalable clustering
- Integrated new binning and MAG evaluation tools: MAGScoT, Galah, Vamb, and MetaDecoder