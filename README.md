# CovidION

Pipeline for analysis of any type or viral organism with small genomes. Pipeline specialized in SARS-CoV-2 MinION samples

For the development of this pipeline it is necessary to download [Guppy](https://community.nanoporetech.com/downloads/guppy/release_notes) from [Nanopore](https://community.nanoporetech.com/downloads) for the basecalling and de-multiplexing steps. If you have used [MinKNOW](https://community.nanoporetech.com/downloads/minion_release/release_notes) to obtain the .fast5 files, you can also perform these steps within the same software.


## Installation of dependencies

If the installation of the yaml file for the envirenment does not correctly install the dependent packages with pip, you can do the installation manually using the following commands:

- [Pangolin](https://github.com/cov-lineages/pangolin): pip install git+https://github.com/cov-lineages/pangolin
- [PangoLEARN](https://github.com/cov-lineages/pangoLEARN.git): pip install git+https://github.com/cov-lineages/pangoLEARN.git
- [Scorpio](https://github.com/cov-lineages/scorpio.git): pip install git+https://github.com/cov-lineages/scorpio.git
- [Constellations](https://github.com/cov-lineages/constellations.git): pip install git+https://github.com/cov-lineages/constellations.git
- [Pango-Designation](https://github.com/cov-lineages/pango-designation.git): pip install git+https://github.com/cov-lineages/pango-designation.git
- [NanoPlot](https://github.com/wdecoster/NanoPlot): pip install NanoPlot && pip install NanoPlot --upgrade
