# Repository to reproduce published analyses of somatic mutations in normal livers 

## Introduction
This repository contains code that was used to generate analyses of somatic mutations in normal livers. A manuscript titled "Somatic mutations and clonal dynamics in healthy and cirrhotic human liver" is currently under review.

## Data availability
Data accompanying this repository is available via Synapse. Many analyses and examples in this repository require the download of the data into folder `liver-pub-repo/data`. The data is available here: `to do: paste public link`.

Please note that we attempted to make as much of the data publicly available as possible. In the case of very large files, the Synapse project contains processed rather than raw data. Raw data is available from EGA (`to do: post repository link here`). Wherever possible we deposited code that should allow to go from raw to processed data. In some cases where intermediate data are large, we deposited a small subset of the data to execute the analysis for that subset. Final variant calls and similar processed data is consistently made available.

## Folder structure:
* lib => shared functions
* analyses => scripts required to reproduce analyses presented in publication