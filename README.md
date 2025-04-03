# Humanatee

A tool for annotating and prioritizing human HiFI-WGS variation.

## Documentation

* [Installation instructions](docs/install.md)
* [Tutorial for TRGT annotation](docs/trgt_tutorial.md)
* [Interpreting TRGT annotation](docs/trgt_interpretation.md)

## Need help?

If you notice any missing features, bugs, or need assistance with analyzing the output of Humanatee,
please don't hesitate to open a GitHub issue.

## Support information

Humanatee is a pre-release software intended for research use only and not for use in diagnostic procedures. While efforts have been made to ensure that Humanatee lives up to the quality that PacBio strives for, we make no warranty regarding this software.
As Humanatee is not covered by any service level agreement or the like, please do not contact a PacBio Field Applications Scientists or PacBio Customer Service for assistance with any Humanatee release. Please report all issues through GitHub instead. We make no warranty that any such issue will be addressed, to any extent or within any time frame.

## Changelog

* 0.1.1
  * Fixed build issue where package data wasn't being copied from source.
  * Fixed bugs that required trgt optional params `--phenotype` and `--gene-lookup` to be specified.
  * Fixed an error in tandem repeat histograms where analysis would stop when the population distribution of alleles is `NoneType`.
  * Fixed error plotting if parental genotypes are missing.

## DISCLAIMER

THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.
