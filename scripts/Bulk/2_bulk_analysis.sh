#!/bin/bash

# Define absolute paths for input and output folders
base_dir="YOUR_PATH/"  # Change this to the absolute path of bulk_analysis
spacer=".."  # Default spacer if needed
jobs_dir="/YOUR_PATH"  # Directory to save .sh job scripts
output_logs="YOUR_PATH/output_logs"  # Directory for output logs
error_logs="YOUR_PATH/error_logs"  # Directory for error logs

# Create jobs_gRNAs directory if it doesn't exist
mkdir -p "$jobs_dir" "$output_logs" "$error_logs"

# Iterate over each folder (g1, g2, ..., g9)
for folder in g1 g2 g3 g4 g5 g6 g9; do
  # Create output_alleles/ directory within each folder (g1, g2, ...)
  output_dir="$base_dir/$folder/output_alleles"
  mkdir -p "$output_dir"
  
  for file in "$base_dir/$folder"/*.csv; do
    filename=$(basename "$file")
    job_name="$jobs_dir/${filename%.*}_alleles_job.sh"  # Save .sh job script in jobs_gRNAs
    
    # Define gRNA reference according to the folder
    if [[ "$folder" == "g1" ]]; then
      myref="GTCCCCTCCACCCCACAGTG"
    elif [[ "$folder" == "g2" ]]; then
      myref="GTACACCCTCAAGCAGTGTG"
    elif [[ "$folder" == "g3" ]]; then
      myref="GTATGCGGATGCAATCTCCG"
    elif [[ "$folder" == "g4" ]]; then
      myref="GTCCCCTCCACCCCACTCTG"
    elif [[ "$folder" == "g5" ]]; then
      myref="CCTCACCTTCCTCACACCGC"
    elif [[ "$folder" == "g6" ]]; then
      myref="GTAGAGTCCCCTCCACCCCACAGTG"
    elif [[ "$folder" == "g9" ]]; then
      myref="ATCCACCTCACCTTCCTCACACCGC"
    else
      echo "No corresponding gRNA found for $folder"
      continue
    fi

    # Create a job script (.sh) for each file and save it in jobs_gRNAs
    cat <<EOT > "$job_name"
#!/bin/bash
#SBATCH --job-name=process_${filename%.*}
#SBATCH --output=$output_logs/${filename%.*}_output.txt  # Save output file in the logs directory
#SBATCH --error=$error_logs/${filename%.*}_error.txt    # Save error file in the logs directory

# Load R module if necessary
#module load R

# Execute the R script with the gRNA as an argument
Rscript /Bulk_3_analysis.R "$file" "$output_dir" "$myref"
EOT

    # Submit the job to SLURM
    sbatch "$job_name"
  done
done
