# Nextflow

Exploring a few ways to parallize this pipeline.

```
jenchang:nextflow $ nextflow run main.nf 
N E X T F L O W  ~  version 21.04.0
Launching `main.nf` [loving_dubinsky] - revision: 405dd752ce
executor >  local (16)
[90/8d8f22] process > mafft (2)    [100%] 8 of 8 ✔
[cc/7f8659] process > fasttree (8) [100%] 8 of 8 ✔
Completed at: 07-Jun-2021 16:59:47
Duration    : 1m 55s
CPU hours   : 0.1
Succeeded   : 16
```

Output directory:

```
Results/
  |_ 01_Align/
  |  |_ PB2_annot_aln.fna
  |  |_ PB1_annot_aln.fna
  |  |_ HH_annot_aln.fna
  |  |_ NP_annot_aln.fna
  |  |_ NS_annot_aln.fna
  |  |_ M_annot_aln.fna
  |  |_ NN_annot_aln.fna
  |  |_ PA_annot_aln.fna
  |
  |_ 02_Trees/
     |_ NN_annot_aln.tre
     |_ PA_annot_aln.tre
     |_ M_annot_aln.tre
     |_ NS_annot_aln.tre
     |_ NP_annot_aln.tre
     |_ PB1_annot_aln.tre
     |_ HH_annot_aln.tre
     |_ PB2_annot_aln.tre
```

To do:

* [ ] automatically draw the trees, color by clade
* [ ] color reassortants in red
