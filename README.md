# ClassificationByNumber
A tool for comparing a nucleotides sequence query  with a database of sequences and identify the closest neighbours that resemble the query.
The CBN uses signal processing and time series techniques to identify similarities between the query sequence and the sequence database.

## Installation
Run the appropriate GCC or Intel installation script provided for each system and the ClassificationByNumbers executable will be generated in folder ./Exec.
```
i.e for a Linux  or a Mac executable with GCC compiler 
First step make the script executable using
chmod u+x ./GCC_Installation

  and then invoke the executable

./GCC_Installation
```
##

## Built With
* [Complete-Striped-Smith-Waterman-Library] (https://github.com/mengyao/complete-striped-smith-waterman-library)
A Smith Waterman alingment algorithm implementation
##

## Usage
```
./ClassificationByNumbers -fa ./References.fasta -r ./ShortReads.fastq -o ./OutputDirectory
```
For more detailed instructions use
```
./ClassificationByNumbers -h
```
## Authors

* Avraam Tapinos
##

## Fundings

This work has been supported by the BBSRC [BB/M001121/1], and the VIROGENESIS project which receives funding from the European Union's Horizon 2020 research and innovation programme under grant agreement No 634650.
##
