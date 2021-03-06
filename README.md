# Python_De_Bruijin_Hw4

## Description
De Bruijin graph assembler operates via prefix sequence to assemble nucleotide strings over a linear (or at least neraly linear) time. Here we implement custom De Bruijin assembler that can be executed from command line.

## Code structure
Basically, code comprises of two parts:
+ code body: object-oriented code comprising of three class objects for graph assembly
+ *argparse* part: integration into command line and code execution

### Code body
The body holds three classes essential for graph assembly
+ *Class **Vertex***: a class for graph vertices which stand for individual k-mers. It contains information about coincindent edges and k-mer coverage
+ *Class **Edge***: a class for graph edges representing existing reads. This class contains information about edge coverage calculated as mean between two interlinked k-mers.
+ *Class **Graph***: the largest classes representing the graph itself. It contains three functions: one for read addition (*add_read*), one for vertex coverage calculation (*calc_init_coverage*) and one for graph visualization (*visualize*) producing a pdf object either with full text (if **t=full**) or only with numeric stats (in any other case).

### Argparse
Arguments for command line go as following:
+ **-i**: a path to fasta file containing reads for further assembly;
+ **-k**: a desired k (k-mer length); default is 3;
+ **-t**: if followed by 'full', draws a full graph with k-mers on vertices and read sequences on edges; if followed by any other random text, produces a stat graph with node and edge coverage
+ **-f** and **-r**: two alternative flags for assembly mode. If **-f** is stated, graph is assembled fron raw read sequences; if **-r** is chosen, graph is assemble from read reverse complements;
+ **function**: a none-flag variable for **t**/**r** alteration handling.

## Launch example
The following command executes graph assembly from forward reads with k equal to 3:
```{bash}
python De_Brujin_Malovichko.py -f -t full -i hw3_dataset.fasta
```
Alternatively, this one will create a stat graph (by random text after **-t**) from reverse reads with k=10:
 ```{bash}
python De_Brujin_Malovichko.py -r -t pop -k 10 -i hw3_dataset.fasta
```

## Output example
Three graphs are stored in this repository as examples for graph assembly: one for forward read assembly in full mode, one for forward read assembly in stat graph mode and the last for full graph assembled from reverse reads. All three were built with default k from reads stored in attached *fasta* file.

## Acknowledgements
Eugene Bakin from Bioinformatics Institute for code backbone and testing material
Pre-existing assemblers for graph visualization function logic
