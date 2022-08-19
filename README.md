# PVC-and-EU-coatomers

This repository contains the structural protein models of the eucaryotic and PVC coatomers (full length and domains, available in ./database), including the comparisons performed with MOMA2 (./MDS_GMM_analyses/input). 

In addition, we performed a clustering analysis of the different types of PVC coatomers to obtain a representative group to compare with the eukaryotic coatomer proteins (./PVC_clustering).

To explore the remote structural relationships present in these proteins, we run a multidimensional scaling analysis and attempt to cluster the results into groups through Gaussian mixture models (./MDS_GMM_analyses/scripts).

Finally, we characterize the structural relationship of a PVC protein (Unicode id A0A225E1Q9) with the EU coatomers (/PVC_example). Firstly, we obtain the structural alignment of this protein with the closest hits. Then, we get a list of orthologous sequences for each hit performed query searches against a Uniprot50 and SwissProt database with jackhmmer (hmmer suite). After that, we obtain their multiple sequence alignment (MSA) with ClustalW to construct sequence blocks, which are aligned in a second step with Mafft using the structural alignment as seed. The resulting MSA show conserved positions that were aligned in the structural comparisons. These MSA can render in a pdf file using ESpript 3.0.
