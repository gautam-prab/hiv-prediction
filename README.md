# hiv-prediction
Predicting Evolution of Resistance Mutations to NRTI Drugs in HIV-1 Using Computational and Statistical Methods.
A science project submitted to the Synopsys Science Fair in 2016, gaining a 2nd prize in biology.

See my paper for a full explanation: https://github.com/gautam-prab/hiv-prediction/blob/master/Final%20Paper.pdf

### Goals
The goal of this project is to create a fitness model for HIV under drug pressure in order to determine which HIV mutations will most likely result from drug stress. The model predicts 22 single mutations that confer resistance to NRTIs and 4,291 pairs of mutations that are not currently in the HIV Resistance Database (http://hivdb.stanford.edu).

### Data
All data was taken from the HIV resistance database in both treated and untreated patients. This study focuses on NRTI drugs, so the Reverse Transcriptase protein was examined.

The datafiles (hiv-untreated.fasta and hiv-treated.fasta) were cut and aligned using UNIPRO UGene into the files (hiv-untreated_subalign.fasta and hiv-treated_subalign.fasta).

The results of gradient descent, resulting in h vectors and J vectors, are stored in hiv-untreated-double-38.hiv and hiv-treated-double-38.hiv, both of which are saved as Java objects. The program compares the sequences to the consensus found at consensus.txt

hiv-subtypes-final-alignment contains alignments of different HIV subtypes but it went unused in this experiment.

The candidate mutations generated by the model are located in FinalData/CandidateSingleMutations.txt and FinalData/CandidateDoubleMutations.txt
