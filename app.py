import streamlit as st
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, molecular_weight

# 1. Configure the Web Layout
st.title("🧬 Molecular Biology Sequence Analyzer")
st.write("An automated tool for quick sequence translation and metrics computation.")

# 2. Text Input Area for Users
dna_input = st.text_area("Enter DNA Sequence (5' to 3')", "ATGCGATCGATCGATCGATCGATC")

# Clean the input sequence data (strip spaces, remove newlines, force uppercase)
clean_dna = dna_input.strip().replace("\n", "").replace(" ", "").upper()

# 3. Validation and Calculation Pipeline
if clean_dna:
    # Define acceptable characters
    valid_nucleotides = set("ATCGN")
    
    # Check if the input contains invalid characters
    if not set(clean_dna).issubset(valid_nucleotides):
        st.error("Invalid characters detected. Please input only A, T, C, G, or N.")
    else:
        # Instantiate a Biopython Seq object
        my_seq = Seq(clean_dna)
        
        st.subheader("📊 Sequence Metrics")
        
        # Calculate GC Content percentage
        gc_content = gc_fraction(my_seq) * 100
        st.metric(label="GC Content", value=f"{gc_content:.2f}%")
        
        st.subheader("🔄 Transcription & Translation")
        
        # Execute transcription and translation
        mrna = my_seq.transcribe()
        protein = my_seq.translate()
        
        st.text_area("mRNA Sequence", str(mrna), height=80)
        st.text_area("Translated Protein Sequence", str(protein), height=80)
        
        st.subheader("🧪 Protein Characteristics")
        try:
            # Remove stop codons (*) prior to molecular weight calculation
            clean_protein = str(protein).replace("*", "")
            if clean_protein:
                prot_weight = molecular_weight(clean_protein, seq_type="protein")
                st.write(f"**Protein Molecular Weight:** {prot_weight:.2f} Da")
            else:
                st.write("**Protein Molecular Weight:** N/A (Empty Sequence)")
        except Exception as e:
            st.warning(f"Could not compute molecular weight: {e}")
