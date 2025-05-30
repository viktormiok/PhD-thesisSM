![](https://img.shields.io/badge/language-R-orange.svg) ![version](https://img.shields.io/badge/GiHub_version-1.1.0-519dd9) ![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/viktormiok/PhD-thesisSM) ![GitHub issues](https://img.shields.io/github/issues/viktormiok/PhD-thesisSM)

![dependencies](https://img.shields.io/badge/dependencies-up%20to%20date-orange)  	![commit](https://img.shields.io/github/last-commit/viktormiok/PhD-thesisSM) ![GitHub](https://img.shields.io/github/license/viktormiok/PhD-thesisSM)

[![Edit with Gitpod](https://gitpod.io/button/open-in-gitpod.svg)](https://gitpod.io/#https://github.com/viktormiok/PhD-thesisSM) 

PhD Thesis Supplementary Materials
========================

This repository includes supplementary materials related to the PhD thesis "Comprehensive molecular characterisation of HPV-induced transformation by longitudinal statistical modelling" defended at VU University Amsterdam.

- [Author](#author)
- [Title](#title)
- [How to compile this thesis](#how-to-compile-this-thesis) 
- [Where to access thesis and supplementary materials](#where-to-access-thesis-and-supplementary-materials)
- [Doctoral thesis supplementary materials](#doctoral-thesis-supplementary-materials)
- [License](#license)

## Author
Viktorian Miok

## Title

**`Comprehensive molecular characterisation of HPV-induced transformation by longitudinal statistical modelling`**

![image](https://user-images.githubusercontent.com/22052679/150094286-6c24c95a-4b20-4807-a269-d63e322be8e2.png)

## How to compile this thesis

Clone this thesis into a local folder using:

```{bash}
clone https://github.com/viktormok/PhD-thesisSM.git
```

To compile the thesis, change to the Thesis folder and use `latexmk`:

```{bash}
cd Thesis
latexmk -pdf thesis.tex
```

## Where to access the thesis and supplementary materials

The pdf version of this thesis can be accessed [here](https://research.vu.nl/en/publications/comprehensive-molecular-characterisation-of-hpv-induced-transform), while the supplementary materials from the publications can be found in this repository.

## Doctoral thesis supplementary materials

I have done the work presented in this thesis in collaboration with the persons listed in each section. 
When using the information in this thesis, please cite the associated papers:


[**Chapter 1:**](https://github.com/viktormiok/PhD-thesisSM/tree/main/Chapter1) `Introduction`

Background, data generation, and statistical analysis.

[**Chapter 2:**](https://github.com/viktormiok/PhD-thesisSM/tree/main/Chapter2) `tigaR: integrative significance analysis of temporal differential gene expression induced by genomic abnormalities`

**Miok, V.**, Wilting, S. M., van de Wiel, M. A., Jaspers, A., van Noort, P. I., Brakenhoff, R. H., Snijders, P. J. F., Steenbergen, R. D. M. and van Wieringen, W. N. [tigaR: integrative significance analysis of temporal differential gene expression induced by genomic abnormalities](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-327), _BMC Bioinformatics_, (2014), 15: 327.

[**Chapter 3:**](https://github.com/viktormiok/PhD-thesisSM/tree/main/Chapter3) `Ridge estimation of the VAR(1) model and its time series chain graph from multivariate time-course omics data`

**Miok, V.**, Wilting, S. M. and van Wieringen, W. N. [Ridge estimation of the VAR(1) model and its time series chain graph from multivariate time-course omics data](https://onlinelibrary.wiley.com/doi/full/10.1002/bimj.201500269), _Biometrical Journal_,(2017), 59(1), 172-191

[**Chapter 4:**](https://github.com/viktormiok/PhD-thesisSM/tree/main/Chapter4) `Ridge estimation of network models from time-course omics data`

**Miok, V.**, Wilting, S.M., van Wieringen, W.N. [Ridge estimation of network models from time‐course omics data](https://onlinelibrary.wiley.com/doi/full/10.1002/bimj.201700195), _Biometrical Journal_, (2019), 61(2), 391-405.

[**Chapter 5:**](https://github.com/viktormiok/PhD-thesisSM/tree/main/Chapter5) `Comprehensive molecular profiling of HPV-induced transformation over time`

Babion, I., **Miok, V.**, Jaspers, A., Huseinovic, A., Steenbergen, R.D., van Wieringen, W.N., Wilting, S.M., [Identification of Deregulated Pathways, Key Regulators, and Novel miRNA-mRNA Interactions in HPV-Mediated Transformation](https://www.mdpi.com/2072-6694/12/3/700), _Cancers_, (2020), 12(3), 700.

[**Chapter 6:**](https://github.com/viktormiok/PhD-thesisSM/tree/main/Chapter6) `Aberrant methylation-mediated silencing of microRNAs contributes to HPV-induced anchorage independence`

Wilting, S.M., **Miok, V.**, Jaspers, A., Boon, D., Sorgard, H., Lando, M., Snoek, B.C., van Wieringen, W.N., Meijer, C.J., Lyng, H. and Snijders, P.J. [Aberrant methylation-mediated silencing of microRNAs contributes to HPV-induced anchorage independence](https://www.oncotarget.com/article/9698/text/), _Oncotarget_, (2016), 7(28), 43805-43819

[**Chapter 7:**](https://research.vu.nl/ws/portalfiles/portal/61554211/chapter+7.pdf) `Conclusion`

Discussion and concluding remarks.

The figures in the introduction are my work. 
When using or adapting them, please cite this thesis as follows:

Viktorian, M. (2019). [_Comprehensive molecular characterisation of HPV-induced transformation by longitudinal statistical modelling_](https://research.vu.nl/en/publications/comprehensive-molecular-characterisation-of-hpv-induced-transform)

## License

__`PhD-thesisSM`__ is distributed under the GPL-3.0 License. Please read the license before using __`PhD-thesisSM`__, which is distributed in the `LICENSE` file.

