
![](logo.png)

![Read the Docs](https://img.shields.io/readthedocs/qmzyme)
[![CircleCI](https://dl.circleci.com/status-badge/img/circleci/LRUEotncYnASivTi54FASD/DM24Fo4B8Af3VgC59vNCCx/tree/main.svg?style=shield)](https://dl.circleci.com/status-badge/redirect/circleci/LRUEotncYnASivTi54FASD/DM24Fo4B8Af3VgC59vNCCx/tree/main)
[![codecov](https://codecov.io/gh/hklem/QMzyme/graph/badge.svg?token=5PISDUT85W)](https://codecov.io/gh/hklem/QMzyme)
[![Powered by MDAnalysis](https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA)](https://www.mdanalysis.org)

[comment]:  <[![PyPI version](https://badge.fury.io/py/QMzyme.svg)](https://badge.fury.io/py/QMzyme)>


*QMzyme is currently under-development. Please note the user interface may change! The first stable API version will be released as QMzyme==1.0.0 on PyPi.*

QMzyme is a Python toolkit to facilitate (quantum mechanical) QM-based enzyme calculations. The GenerateModel module guides the process of generating calculation ready truncated or partitioned molecule regions. Any input file(s) accepted by [MDAnalysis to create a Universe](https://userguide.mdanalysis.org/stable/universe.html) object can be used to start working in QMzyme. From there, the code relies on more flexible QMzyme objects: QMzymeAtom, QMzymeResidue, QMzymeRegion and QMzymeModel. 

Full documentation with installation instructions, technical details and examples can be found in [Read the Docs](https://qmzyme.readthedocs.io/).

## Contributing to QMzyme
For suggestions and improvements of the code (greatly appreciated!), please reach out through the issues and pull requests options of Github.  

## References

QMzyme has not been formally published. Please refer back here for the proper citation once it has.

QMzyme's main dependency is MDAnalysis. When using QMzyme please [cite MDAnalysis accordingly](https://www.mdanalysis.org/pages/citations/):

* R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler, D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein. MDAnalysis: A Python package for the rapid analysis of molecular dynamics simulations. In S. Benthall and S. Rostrup, editors, Proceedings of the 15th Python in Science Conference, pages 98-105, Austin, TX, 2016. SciPy, doi:10.25080/majora-629e541a-00e.

* N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein. MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319-2327, doi:10.1002/jcc.21787. PMCID:PMC3144279

If you use QMzyme to write QM software calculation input files, please include this citation for [AQME](https://aqme.readthedocs.io/en/latest/):  

* Alegre-Requena, J. V.; Sowndarya, S.; Pérez-Soto, R.; Alturaifi, T.; Paton, R. AQME: Automated Quantum Mechanical Environments for Researchers and Educators. Wiley Interdiscip. Rev. Comput. Mol. Sci. 2023, 13, e1663. (DOI: 10.1002/wcms.1663).  


#### Copyright
Copyright (c) 2024, Heidi Klem

#### Acknowledgements
Project architecture initially inspired by [MolSSI CMS Cookiecutter](https://github.com/molssi/cookiecutter-cms).
