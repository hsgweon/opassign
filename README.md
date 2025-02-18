# OpAssign: Taxonomic assigner on 16S-ITS-23S operons using RDP Classifier.

Taxonomic assignment on long read sequences (ONT/PacBio) using RDP Classifier trained from GROND database.

---

## ✨ Features

✅ Assign taxonomy

✅ Uses RDP Classifier to classify. 

✅ Supports multiple CPUs

## ⚡ Installation

### Clone the Repository & Set Up the Environment

```bash
# Clone the repo
git clone https://github.com/hsgweon/opassign.git

# Create the seqdemu environment (ensure that you have conda installed)
mamba create -n opassign_env -y -c conda-forge -c bioconda conda-forge::biopython conda-forge::psutil progressbar2 
```

## 🎬 Running seqdemu

```bash
mamba activate opassign_env
opassign.py -i rep_16S23S_nr.fasta -o assigned_tanoxony.txt -t 100
```

