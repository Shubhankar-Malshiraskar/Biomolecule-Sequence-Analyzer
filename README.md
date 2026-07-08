# Biomolecule Sequence Analyser 

A fully deployable, modular bioinformatics dashboard built with Python. This tool bridges the gap between raw sequence data and actionable laboratory metrics, providing an automated pipeline for sequence validation, transcription, translation, and downstream analysis.

Designed with a focus on clean user experience (UX) and data chunking, the interface mitigates cognitive fatigue common in traditional, clustered bioinformatics software.

## Biological Utility
This application automates routine calculations required for wet-lab preparation and microbiological analysis:
* **Assay Preparation:** Calculates exact molecular weights of translated peptide chains, crucial for standardizing protein extraction quantifications (e.g., Lowry assays).
* **PCR Optimization:** Computes highly precise Nearest Neighbor (NN) melting temperatures ($T_m$) and GC content for custom primers, bypassing manual Wallace rule approximations.
* **Data Integrity:** Features strict input validation and sanitization for DNA, RNA (including mRNA and rRNA), and standard amino acid sequences.

## UI Architecture & Tech Stack
The application is structured as a single-page dashboard utilizing modular tabs to separate concerns, ensuring that complex data outputs (like alignments) do not clutter primary metrics. 

**Dependencies:**
* **Streamlit:** Handles the frontend UI, session state management, and real-time cloud deployment.
* **Biopython:** Powers the core biological logic (`Bio.Seq`, `Bio.SeqUtils.MeltingTemp`, `Bio.Align`).
* **Pandas:** Manages internal data structures and CSV export formatting.

## Key Features
* **Session State Memory:** Automatically tracks and logs analyzed sequences during an active session, allowing users to export their query history as a `.csv` file.
* **Dynamic Visualizations:** Generates immediate nucleotide/amino acid frequency distribution charts based on the sanitized input.
* **Global Sequence Alignment:** Implements the Needleman-Wunsch algorithm via `Bio.Align` for pairwise target matching.
* **Target Binding Check:** Automatically scans the primary sequence and its reverse complement to verify primer binding loci.

## Local Deployment (For Developers)
If you wish to run this tool locally bypassing the Streamlit Community Cloud:

1. Clone this repository:
   ```bash
   git clone [https://github.com/your-username/sequence-analyzer.git](https://github.com/your-username/sequence-analyzer.git)
