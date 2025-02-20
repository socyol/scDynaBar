#!/bin/bash

#SBATCH --job-name=blast_job          
#SBATCH --output=blast_output_%j.txt   
#SBATCH --error=blast_error_%j.txt     


#######      SLURM JOB         ###########
##########################################


# Fasta file with sequences to analyze with BLAST
fasta_file="scDynaBar/scripts/SingleCells/Blast_function/gRNA_sequences_with_frequencies_DAY4-10.fasta"

# Output file
output_file="scDynaBar/scripts/SingleCells/Blast_function/blast_results_DAY4_10.txt"

# BLAST database
blast_db="mm10_db"

# Clean outputfile in case it exists
> "$output_file"



sequence=""


while read -r line; do
    if [[ $line == ">"* ]]; then
       
        if [[ -n "$sequence" ]]; then
            echo "---------------------" >> "$output_file"
            echo "$sequence" > temp_query.fasta
            blastn -query temp_query.fasta -db "$blast_db" -outfmt 6 >> "$output_file"
            echo "---------------------" >> "$output_file"
            sequence=""  
        fi
        
        echo "$line" >> "$output_file"
    else
        
        sequence+="$line"
    fi
done < "$fasta_file"



if [[ -n "$sequence" ]]; then
    echo "---------------------" >> "$output_file"
    echo "$sequence" > temp_query.fasta
    blastn -query temp_query.fasta -db "$blast_db" -outfmt 6 >> "$output_file"
    echo "---------------------" >> "$output_file"
fi


rm temp_query.fasta

echo "Results of BLAST in $output_file"

