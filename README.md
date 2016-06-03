# STM-seq-count
Perl scripts for STM-NGS sequence tag counting.

The current project has two perl scripts to be used as a pipeline to count for ocurrencies  of sets of transposon signature tags within next generation sequencing (NGS) reads. 

These scripts are written for STM experiments based on Sinorhizobium meliloti signature-tagged mini-Tn5 transposon library, but can be adapted for other cases.

Pobigaylo, N., Wetter, D., Szymczak, S., Schiller, U., Kurtz, S., Meyer, F., Nattkemper, T.W., Becker, A., 2006. Construction of a large signature-tagged mini-Tn5 transposon library and its application to mutagenesis of Sinorhizobium meliloti. Appl Env. Microbiol 72, 4329–4337.

The first script, count_STM.pl, works by loading in the memory each sequence read, in a one by one fashion, and then looks for the sequence signatures corresponding to test condition/replicates, pool, KpnI or HindIII restriction site, and mutant-tag; generating a text file containing on each line the count of any condition/pool/mutant-tag found. 

The second script, order_STM.pl, searches for replicates, input and output conditions and generates two text files (corresponding to input and output conditions respectively) containing the count of all the replicates for each condition/pool/mutant-tag. 

Next (not included on this project), a standard SQL database software has to be used to join input and output conditions. As a final step, normalization and differential expression analysis has to be done by third party software. 
