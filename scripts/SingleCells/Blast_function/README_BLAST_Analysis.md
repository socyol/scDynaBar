# 🔬 Off-Target Analysis of gRNAs with BLAST

This repository contains scripts to analyze **off-target effects of gRNAs** using BLAST in a CRISPR barcoding experiment (Day 4 and Day 10). The analysis consists of two main steps:

1. **`1_blast_analysis.R`**: Generates a FASTA file with gRNA sequences and a CSV file with their frequencies.
2. **`2_blast_analysis.sh`**: Runs BLAST on a computational cluster to identify potential off-target sites.

---

## 📂 Project Structure

```
/Blast_function
│── 1_blast_analysis.R          # R script to generate FASTA and CSV files
│── 2_blast_analysis.sh         # SLURM script to execute BLAST analysis
│── barcodes-alleles_cuts_day4-10.csv  # Input: gRNA sequences and counts
│── gRNA_sequences_with_frequencies_DAY4-10.fasta  # Output: FASTA file for BLAST
│── gRNA_sequences_with_frequencies_DAY4-10.csv  # Output: CSV with sequence frequencies
│── blast_results_DAY4_10.txt   # Final BLAST output
```

---

## 🚀 **Step 1: Generating FASTA & CSV Files with gRNAs**
📌 **Script:** `1_blast_analysis.R`

This script reads a CSV file containing **gRNA sequences**, calculates their **relative frequency**, and generates:

- A **FASTA file** (`.fasta`) for BLAST.
- A **CSV file** (`.csv`) with gRNA sequences and their frequency.

### **🔧 Usage**
Run the following command in R:

```r
source("1_blast_analysis.R")
```

### **📌 Input**
- `barcodes-alleles_cuts_day4-10.csv`: Contains **gRNA sequences** from the experiment.

### **📌 Output**
- `gRNA_sequences_with_frequencies_DAY4-10.fasta`: FASTA file for BLAST.
- `gRNA_sequences_with_frequencies_DAY4-10.csv`: CSV file with gRNA frequencies.

---

## 🚀 **Step 2: Running BLAST on a Cluster**
📌 **Script:** `2_blast_analysis.sh`

This SLURM script processes the **FASTA file** using *BLAST* to search for **off-target sites** in the **mouse genome (mm10)**.

### **🔧 Submitting the Job**
Run the following command in the cluster:

```bash
sbatch 2_blast_analysis.sh
```

### **📌 Input**
- `gRNA_sequences_with_frequencies_DAY4-10.fasta` (Generated from Step 1)
- `mm10_db` (Mouse genome BLAST database)

### **📌 Output**
- `blast_results_DAY4_10.txt`: BLAST results in tabular format.

---

## 🔎 **BLAST Execution Details**
This script processes the **FASTA file** line by line, running **BLASTN** for each sequence. It:
1. Reads the **FASTA file** (`.fasta`).
2. Extracts **gRNA sequences**.
3. Runs **BLASTN** against the **mouse genome (`mm10_db`)**.
4. Saves results in `blast_results_DAY4_10.txt`.

### **Example of a BLAST result (tabular format):**
```
gRNA1	chr1	98.7	23	0	0	100	123	+	1e-05	45.6
gRNA2	chr3	95.2	23	1	0	210	233	-	1e-03	42.1
```
Columns:
1. gRNA ID
2. Chromosome
3. Identity %
4. Alignment length
5. Mismatches
6. Gaps
7. Start position
8. End position
9. Strand (+/-)
10. E-value
11. Bit score

---

## 📌 **Transferring Files to the Cluster**
To transfer files to the cluster, use `scp`:

```bash
scp gRNA_sequences_with_frequencies_DAY4-10.fasta user@server:/path/to/cluster/
scp 2_blast_analysis.sh user@server:/path/to/cluster/
```

To copy results back:

```bash
scp user@server:/path/to/cluster/blast_results_DAY4_10.txt .
```

---

## 📢 **Final Notes**
- The **FASTA file** must be correctly formatted before running BLAST.
- The **BLAST database (`mm10_db`)** should be indexed in the cluster.
- Check **BLAST logs (`blast_output_*.txt` and `blast_error_*.txt`)** for any errors.

🚀 **For any questions, feel free to reach out!** ✨
