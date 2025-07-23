PTTseek is a python script designed to detect change points in RNA-seq data and has been successfully applied to metabolically labelled RNA-seq (e.g. 4sU-seq) for detecting premature termination sites. It uses bigwig and bed file as input and the output is a tsv file with genomic coordinates of drop sites and the corresponding drop values.

Requirements:

python3
numpy
pyBigWig
ruptures
argparse

Command for running the script:

python3 PTTseek.py -b your_bigwig_file -i your_bed_file -o output_file.tsv --bin 10 --min-length 100
