# Genome Analysis Pipelines

**Decontamination** – A multi-step workflow for detecting possible contamination in genome assemblies based on protein annotations.

This repository contains some important things you should know before running this **contamination detection** script:

## 1. Decontamination

The **Decontamination** pipeline is a **Python + Diamond workflow** to identify possible contamination in genome sequencing data.
It integrates similarity searches against curated databases and KEGG functional assignments.

## Quick Start

1. Prepare your files:
   - `braker.aa` predicted proteome
   - Reference proteome (type species)
   - Diamond databases: SwissProt and NR
   - KEGG GhostKOALA account

2. Run pipeline:

```bash
# Reference comparison
python3 1_dmnd.py \
-f /path/to/predicted_proteome.aa \
-d /path/to/diamond_db.dmnd \
-t 15 \
-o strain_reference

python3 2_dmnd_header.py \
-f /path/to/predicted_proteome.aa \
-d /path/to/_Dblastp.txt \
-o strain

# Extract nohits
python3 3_extraction.py \
-b /path/to/header_outfmt7.txt \
-p /path/to/predicted_proteome.aa \
-o strain


# Align nohits against SwissProt and NR
python 1_dmnd.py nohits.faa swissprot.dmnd
python 2_dmnd_header.py
python 1_dmnd.py nohits.faa nr.dmnd
python 2_dmnd_header.py

# KEGG GhostKOALA
# (Upload nohits.faa online [https://www.kegg.jp/ghostkoala/], download .top result)

# Score KEGG
python 4_score_kegg.py

# Combine results
python3 5_intersection.py \
-s /path/to/strain_sprot_outfmt7.txt \
-n /path/to/strain_nr_outfmt7.txt \
-k /path/to/strain_kegg.top \
-ks 50 \
-f /path/to/strain_nohits.faa \
-o strain


### Workflow
                ┌──────────────────────┐
                │  Predicted proteome  │
                │ (BRAKER3 .aa output) │
                └─────────┬────────────┘
                          │
                          ▼
              ┌──────────────────────────┐
              │  1_dmnd.py vs reference  │
              └─────────┬────────────────┘
                        │
                        ▼
              ┌──────────────────────────┐
              │     2_dmnd_header.py     │
              └─────────┬────────────────┘
                        │
                        ▼
              ┌──────────────────────────┐
              │ 3_extraction.py          │
              │ → produces nohits.faa    │
              └─────────┬────────────────┘
                        │
        ┌───────────────┼───────────────────────┐
        ▼               ▼                       ▼
┌──────────────────┐  ┌───────────────┐    ┌─────────────────┐
│ 1+2 vs SwissProt │  │   1+2 vs NR   │    │Upload nohits.faa│
│  diamond+sprot   │  │  diamond+nr   │    │  to KEGG KOALA  │
└───────┬──────────┘  └───────┬───────┘    └────────┬────────┘
        ▼                     ▼                     ▼
        └────────────────┬────┴──────────────┬──────┘
                        ▼                    ▼
              ┌───────────────────┐
              │ 4_score_kegg.py   │
              └────────────────┬──┘
                               ▼
              ┌─────────────────────┐
              │  5_intersection.py  │
              │ → intersection.csv  │
              └─────────────────────┘



### Prerequisites

**Files:**

* Predicted proteome (tested with **BRAKER3** output `braker.aa`).
  *Suggestion: filter the genome and remove contigs <500 bp before prediction.*
* Predicted proteome of the **type species** of the genus under study.
* **Updated NR** database indexed for Diamond.
* **Updated SwissProt** database indexed for Diamond.

**Programs:**

* [Diamond](https://github.com/bbuchfink/diamond) (for fast sequence alignment).
* KEGG GhostKOALA (online tool, requires institutional email. https://www.kegg.jp/ghostkoala/).
* Python v3.
* Linux environment.

---

### Workflow

The pipeline is organized in **five Python scripts**. Below is the recommended workflow:

1. **`1_dmnd.py`**
   Run Diamond to align your predicted proteome against a reference proteome.

   ```bash
   python 1_dmnd.py my_genome_braker.aa reference_proteome.faa
   ```

2. **`2_dmnd_header.py`**
   Parse Diamond output to clean and format headers.

3. **`3_extraction.py`**
   Extract unmatched proteins and generate a `nohits.faa` file for further analysis.

4. Run Diamond again on the `nohits.faa` (use the **`1_dmnd.py`**):

   * Against **SwissProt**
   * Against **NR**

   After each run, repeat step 2 (**`2_dmnd_header.py`**) on the outputs.

5. Upload `nohits.faa` to **KEGG GhostKOALA**, download the `.top` result, and place it in your working directory.

6. **`4_score_kegg.py`**
   Combine KEGG results with one of the Diamond outputs to calculate a suggested **KEGG score cutoff**.

7. **`5_intersection.py`**
   Merge NR, SwissProt, and KEGG results into an **intersection table**.

---

### Output

* `nohits.faa` – proteins without hits in reference comparison.
* Formatted Diamond result tables.
* KEGG score cutoff suggestion file.
* Final intersection table (`intersection.csv`) – ready for import into visualization tools.

---

### Manual inspection

The final step is to import the **intersection tables** into your preferred visualization or analysis tool (e.g., R, Excel, Python notebooks, Google Sheets) to manually assess possible contaminations.

---

## Citation

Coming soon.

---

## About the author

**Dáfne de Oliveira Vianei**
Biological Sciences, Universidade Federal de Minas Gerais
📧 [vianeidafne@gmail.com](mailto:vianeidafne@gmail.com)
🔗 [LinkedIn](https://www.linkedin.com/in/dafnevianei)

