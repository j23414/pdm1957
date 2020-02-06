# beast

Any kind of Bayesean analysis appreciated.

```
 beast/
  |_HH_strictsubset.fna      # Had to subset down to only H2, because not enough data to fill in tree 
  |_strict.tre               # Tree
  |
  | # Following beast analysis provided by Michael Z.
  |_pop_size.xml             # feed this to beast, I think this is finding effective population size?
  |_state_transitions.xml    # also feed this to beast, I think this is regional state transitions?

```

## Installation and running on MacOS


```
# Install
brew install beast

# Run a beauti generated xml file
beast pop_size.xml
```

## Description of Models Used

**Beast Primary Citation**

* Suchard MA, Lemey P, Baele G, Ayres DL, Drummond AJ, Rambaut A. [Bayesian phylogenetic and phylodynamic data integration using BEAST 1.10.](https://pubmed.ncbi.nlm.nih.gov/29942656) Virus evolution. 2018 Jan;4(1):vey016.



**Uncorrelated relaxed clock**

* Drummond AJ, Ho SY, Phillips MJ, Rambaut A. [Relaxed phylogenetics and dating with confidence.](https://pubmed.ncbi.nlm.nih.gov/16683862) PLoS biology. 2006 May;4(5).

**GTR nucleotide substitution model**

* Tavare S. [Some probabilistic and statistical problems in the analysis of DNA sequences.](http://www.damtp.cam.ac.uk/user/st321/CV_&_Publications_files/STpapers-pdf/T86.pdf) Lectures on mathematics in the life sciences. 1986 Dec 31;17(2):57-86.

**Discrete gamma-distributed rate heterogeneity model**

* Yang Z. [Maximum likelihood phylogenetic estimation from DNA sequences with variable rates over sites: approximate methods.](https://pubmed.ncbi.nlm.nih.gov/7932792) Journal of Molecular evolution. 1994 Sep 1;39(3):306-14.

**Skyride Coalescent, tree density model**

* Minin VN, Bloomquist EW, Suchard MA. [Smooth skyride through a rough skyline: Bayesian coalescent-based inference of population dynamics.](https://pubmed.ncbi.nlm.nih.gov/18408232) Molecular biology and evolution. 2008 Jul 1;25(7):1459-71.

**CTMC Scale Reference Prior Model**

* Ferreira MA, Suchard MA. [Bayesian analysis of elapsed times in continuous‐time Markov chains.](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.564.4652&rep=rep1&type=pdf) Canadian Journal of Statistics. 2008 Sep;36(3):355-68.



## To Do

* Create a fasta2xml.sh script to loop over all 8 gene segments.
* Create a checkBeastInput.sh script to determine if the tree is too sparse or gappy to coalesce.
