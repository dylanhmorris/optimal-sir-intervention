# Optimal, near-optimal, and robust epidemic control
[Dylan H. Morris](https://dylanhmorris.com)(1\*), [Fernando W. Rossine](https://scholar.princeton.edu/ctarnita/people/fernando-rossine)(1\*), [Joshua B. Plotkin](https://www.bio.upenn.edu/people/joshua-plotkin)(2), and [Simon A. Levin](https://slevin.princeton.edu/)(1).

\* These authors contributed equally

1. [Dept. of Ecology and Evolutionary Biology](http://eeb.princeton.edu/), Princeton University, Princeton, NJ, USA
3. [Dept. of Biology](https://www.bio.upenn.edu/) \& [Dept. of Mathematics](https://www.math.upenn.edu/), The University of Pennsylvania, Philadelphia, PA, USA

## Repository information
This repository accompanies the article "Optimal, near-optimal, and robust epidemic control" (Morris, Rossine et al.). It provides code for reproducing numerical solutions of equations in the paper and for recreating associated display figures.

## License and citation information
If you use the code or data provided here, please make sure to do so in light of the project [license](LICENSE.txt) and please cite our work as below:

- D.H. Morris, F.W. Rossine et al. Optimal, near-optimal, and robust epidemic control. 2020. https://arxiv.org/abs/2004.02209.

Bibtex record:
```
@article{morris2021optimal,
    Author = {Morris, Dylan H. and
    Rossine, Fernando W. and 
    Plotkin, Joshua B. and
    Levin, Simon A.},
    Title = {Optimal, near-optimal, and robust epidemic control},
    Date = {2020},
    journal = {arXiv and OSF preprint},
    doi = {10.31219/osf.io/9gr7q},
    url = {https://arxiv.org/abs/2004.02209}
}
```

## Article abstract 
In the absence of drugs and vaccines, policymakers use non-pharmaceutical interventions such as social distancing to decrease rates of disease-causing contact, with the aim of reducing or delaying the epidemic peak. These measures carry social and economic costs, so societies may be unable to maintain them for more than a short period of time. Intervention policy design often relies on numerical simulations of epidemic models, but comparing policies and assessing their robustness demands clear principles that apply across strategies. Here we derive the theoretically optimal strategy for using a time-limited intervention to reduce the peak prevalence of a novel disease in the classic Susceptible-Infectious-Recovered epidemic model. We show that broad classes of easier-to-implement strategies can perform nearly as well as the theoretically optimal strategy. But neither the optimal strategy nor any of these near-optimal strategies is robust to implementation error: small errors in timing the intervention produce large increases in peak prevalence. Our results reveal fundamental principles of non-pharmaceutical disease control and expose their potential fragility. For robust control, an intervention must be strong, early, and ideally sustained.

## Directories
- ``src``: all code, including numerics and figure generation:
- ``out``: output files
    - ``out/results``: numerical analysis outputs as comma-separated values (``.csv``) files. 
    - ``out/figures``: figures generated from results

## Reproducing analysis

A guide to reproducing the analysis from the paper follows.

### Getting the code
First download the code. The recommended way is to ``git clone`` our Github repository from the command line:

    git clone https://github.com/dylanhmorris/optimal-sir-intervention.git

Downloading it manually via Github's download button or on OSF should also work.

### Dependency installation
The analysis can be auto-run from the project ``Makefile``, but you may need to install some external dependencies first. In the first instance, you'll need a working installation of Python 3 with the package manager ``pip3`` and a working installation of Gnu Make or similar. A few external python packages can then be installed by typing.

    make depend

You may also need a working TeX installation to render the text for the figures. If you do not have TeX, you can get around this by setting ``mpl.rcParams['text.usetex'] = False`` in ``plotting_style.py``.

### Running the analyses

The simplest approach is simply to type ``make`` at the command line, which should produce a full set of figures and results.

If you want to do things piecewise, typing ``make <filename>`` for any of the files present in the complete repository uploaded here should also work.

Some shortcuts are available:

- ``make results`` calculates numerical results that are slower to calculate (and therefore must be saved to disk)
- ``make figures`` produces all figures
- ``make clean`` removes all generated files, leaving only source code (though it does not uninstall packages)
