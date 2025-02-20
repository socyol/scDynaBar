# 📊 Bulk Analysis Pipeline

This repository contains scripts for processing bulk sequencing data from CRISPR experiments. The analysis consists of three main steps:

1. **`1_bulk_analysis.R`**: Processes FASTQ files into CSV format, extracting sequence and quality information.
2. **`2_bulk_analysis.sh`**: Generates SLURM job scripts for processing gRNA sequences and submits them to the cluster.
3. **`3_bulk_analysis.R`**: Executes the final R analysis for each gRNA, using reference sequences.

---

## 📂 Project Structure

```
/Bulk_analysis
│── 1_bulk_analysis.R          # R script to process FASTQ files into CSVs
│── 2_bulk_analysis.sh         # SLURM script to process CSV files and run Bulk_3_analysis.R
│── Bulk_3_analysis.R          # R script for final analysis of gRNA sequences
│── FASTQ_FOLDER/              # Directory containing raw FASTQ files
│── Processed_CSVs/            # Directory where processed CSVs are saved
│── output_logs/               # Directory for SLURM output logs
│── error_logs/                # Directory for SLURM error logs
│── gRNA_Folders/              # Directory containing results grouped by gRNA
```

---

## 🚀 **Step 1: Processing FASTQ Files**
📌 **Script:** `1_bulk_analysis.R`

This script processes raw **FASTQ files** and extracts:
- **Nucleotide sequences**
- **Quality scores**
- **Read identifiers**

### **🔧 Usage**
Modify `fastq_dir` and `output_dir` in the script, then run:
```r
source("1_bulk_analysis.R")
```

### **📌 Input**
- FASTQ files stored in `FASTQ_FOLDER/`

### **📌 Output**
- CSV files with sequence data stored in `Processed_CSVs/`

---

## 🚀 **Step 2: Generating and Running SLURM Jobs**
📌 **Script:** `2_bulk_analysis.sh`

This script:
1. Creates **SLURM job scripts** for processing gRNA sequences.
2. Assigns a **reference gRNA sequence** to each dataset.
3. Submits **jobs to the cluster** using `sbatch`.

### **🔧 Submitting Jobs**
Modify `base_dir` in the script, then run:
```bash
bash 2_bulk_analysis.sh
```

### **📌 Input**
- CSV files from `Processed_CSVs/`

### **📌 Output**
- Processed CSV files in `gRNA_Folders/`
- SLURM logs in `output_logs/` and `error_logs/`

---

## 🚀 **Step 3: Running Bulk_3_analysis.R**
📌 **Script:** `Bulk_3_analysis.R`

Each SLURM job calls this script, which:
- Reads the **input CSV file**.
- Processes **allele sequences**.
- Saves results in **output_alleles/**.

This script runs automatically when `2_bulk_analysis.sh` is executed.

---

## 📌 **File Transfers to the Cluster**
To transfer files to the cluster, use:
```bash
scp -r Processed_CSVs/ user@server:/path/to/cluster/
scp 2_bulk_analysis.sh user@server:/path/to/cluster/
```

To copy results back:
```bash
scp -r user@server:/path/to/cluster/gRNA_Folders/ .
```

---

## 📢 **Final Notes**
- Make sure to edit `YOUR_PATH` in the scripts before running.
- Check SLURM logs in `output_logs/` and `error_logs/` for any issues.
- Modify reference gRNA sequences in `2_bulk_analysis.sh` if needed.

🚀 **For any questions, feel free to reach out!** ✨
