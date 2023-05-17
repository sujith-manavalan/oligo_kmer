# Oligo_Kmer Tool

The rapid detection of pathogenic microbes is crucial for biosurveillance, disease diagnostics, and epidemiological studies. Oligonucleotide-based hybridization methods provide cost-effective and highly specific detection even at the strain level.

The Oligo_Kmer tool is designed to build a database of unique oligonucleotides for pathogens, which can be used for designing hybridization chips for the rapid identification of pathogens.

## Features

- **Customization in Oligo Design**: Users can specify their own parameters for oligo design accprding to the chip design needs including oligo size, Hamming distance threshold, and melting temperature range.

- **Handles Large Oligo Sizes**: The tool has been optimized to handle large oligo sizes.

- **SQL Database Storage**: Oligo data is stored in a MySQL database, facilitating global collaborative efforts.

- **GUI Interface**: A GUI interface is provided for easy data retrieval from the SQL database, enhancing broader user accessibility.

## How to Use

Please note: Before running the script, edit the MySQL database parameters in the code to match your database configuration.

To run the `Oligo_Kmer` python script, you will need to provide the following arguments. 

- `-f`: Path to the multi-FASTA file.
- `-k`: Required oligomer size.
- `-d`: Hamming distance threshold.
- `-min_tm`: Minimum Tm in Celsius.
- `-max_tm`: Maximum Tm in Celsius.

A small fasta file(mini.fasta) is provided for testing 

Example useage:

```bash
python oligo_kmer.py -f /path/to/your/file.fasta -k 25 -d 5 -min_tm 50 -max_tm 60
```

## Future Developments

I welcome any suggestions/contributions for improvement of the code . I am hoping to further optimize the code for larger datasets, add multi-threading feature as well as additional customization options.
