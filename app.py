import streamlit as st
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, molecular_weight
from Bio.SeqUtils import MeltingTemp as mt
from Bio import Align
import pandas as pd
import io

# 1. Page Configuration
st.set_page_config(page_title="Biomolecule Sequence Analyser", layout="wide")

# Initialize session state for sequence history
if "sequence_history" not in st.session_state:
    st.session_state.sequence_history = []

# --- SIDEBAR: Settings, Input, and History ---
with st.sidebar:
    st.title("Biomolecule Sequence Analyser")
    
    # Sequence Type Selector
    seq_type = st.selectbox(
        "Sequence Type", 
        ["DNA", "RNA (General)", "mRNA", "rRNA", "Protein"]
    )
    
    # Dynamic placeholder based on selected type
    if seq_type == "DNA":
        default_seq = "ATGCGATCGATCGATCGATCGATC"
    elif "RNA" in seq_type:
        default_seq = "AUGCGAUCGAUCGAUCGAUCGAUC"
    else:
        default_seq = "MKWVTFISLLFLFSSAYSRGVFRR"
    
    # Main Input
    raw_input = st.text_area("Paste your primary sequence here", default_seq, height=150)
    clean_seq = raw_input.strip().replace("\n", "").replace(" ", "").upper()
    
    if st.button("Analyze Sequence"):
        if clean_seq and clean_seq not in st.session_state.sequence_history:
            st.session_state.sequence_history.insert(0, clean_seq)
            st.session_state.sequence_history = st.session_state.sequence_history[:10] # Keep last 10
            
    st.markdown("---")
    
    # History & Export
    st.subheader("Session History")
    if not st.session_state.sequence_history:
        st.caption("No history in this session.")
    else:
        for i, seq in enumerate(st.session_state.sequence_history):
            st.text(f"{i+1}. {seq[:15]}...")
            
        # CSV Export Feature
        history_df = pd.DataFrame(st.session_state.sequence_history, columns=["Sequence"])
        csv = history_df.to_csv(index=False).encode('utf-8')
        st.download_button(
            label="Download History (CSV)",
            data=csv,
            file_name='sequence_history.csv',
            mime='text/csv',
        )

# --- MAIN PAGE: Modular Analysis Tabs ---
if not clean_seq:
    st.info("Please enter a sequence in the sidebar to begin.")
else:
    # Validation Logic
    is_valid = True
    error_msg = ""
    
    if seq_type == "DNA":
        valid_chars = set("ATCGN")
        if not set(clean_seq).issubset(valid_chars):
            is_valid = False
            error_msg = "Invalid characters. Input only A, T, C, G, or N."
    elif "RNA" in seq_type:
        valid_chars = set("AUCGN")
        if not set(clean_seq).issubset(valid_chars):
            is_valid = False
            error_msg = "Invalid characters. Input only A, U, C, G, or N."
    else:
        valid_chars = set("ACDEFGHIKLMNPQRSTVWY*")
        if not set(clean_seq).issubset(valid_chars):
            is_valid = False
            error_msg = "Invalid characters. Input valid standard amino acids."

    if not is_valid:
        st.error(error_msg)
    else:
        my_seq = Seq(clean_seq)
        
        # 5 Clean, Modular Tabs
        tab_overview, tab_visuals, tab_conversions, tab_align, tab_primer = st.tabs([
            "Overview", 
            "Visualizations", 
            "Conversions",
            "Sequence Aligner",
            "Primer Analysis"
        ])
        
        # --- TAB 1: Overview ---
        with tab_overview:
            st.markdown(f"### {seq_type} Core Metrics")
            if seq_type == "DNA" or "RNA" in seq_type:
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric(label="Total Length", value=f"{len(clean_seq):,} bp")
                with col2:
                    gc_content = gc_fraction(my_seq) * 100
                    st.metric(label="GC Content", value=f"{gc_content:.1f}%")
                with col3:
                    st.metric(label="AT/AU Content", value=f"{(100 - gc_content):.1f}%")
            else:
                col1, col2 = st.columns(2)
                with col1:
                    st.metric(label="Total Length", value=f"{len(clean_seq):,} AA")
                with col2:
                    try:
                        clean_protein = clean_seq.replace("*", "")
                        st.metric(label="Estimated Mass", value=f"{molecular_weight(clean_protein, seq_type='protein'):,.2f} Da")
                    except Exception:
                        st.metric(label="Estimated Mass", value="N/A")
            
            st.markdown("---")
            with st.expander("View Raw Cleaned Sequence"):
                st.code(clean_seq, language="text")

        # --- TAB 2: Visualizations ---
        with tab_visuals:
            if seq_type == "DNA":
                counts = {"A": clean_seq.count("A"), "T": clean_seq.count("T"), "C": clean_seq.count("C"), "G": clean_seq.count("G")}
            elif "RNA" in seq_type:
                counts = {"A": clean_seq.count("A"), "U": clean_seq.count("U"), "C": clean_seq.count("C"), "G": clean_seq.count("G")}
            else:
                counts = {aa: clean_seq.count(aa) for aa in set(clean_seq) if aa in valid_chars and aa != "*"}
                counts = dict(sorted(counts.items()))
                
            st.bar_chart(pd.DataFrame(list(counts.values()), index=list(counts.keys()), columns=["Count"]), height=300)

        # --- TAB 3: Conversions ---
        with tab_conversions:
            if seq_type == "DNA":
                col1, col2 = st.columns(2)
                with col1:
                    st.markdown("**Transcription (mRNA)**")
                    st.code(str(my_seq.transcribe()), language="text")
                with col2:
                    st.markdown("**Translation (Protein)**")
                    st.code(str(my_seq.translate()), language="text")
            elif "RNA" in seq_type:
                st.markdown("**Translation (Protein)**")
                st.code(str(my_seq.translate()), language="text")
            else:
                st.write("Sequence is already a protein. Reverse translation is not computationally absolute.")

        # --- TAB 4: Pairwise Sequence Aligner ---
        with tab_align:
            st.markdown("### Pairwise Global Alignment")
            st.write("Compare the primary sequence against a target sequence to find regions of similarity.")
            target_seq_raw = st.text_area("Input Target Sequence for Alignment")
            target_seq = target_seq_raw.strip().replace("\n", "").replace(" ", "").upper()
            
            if st.button("Run Alignment"):
                if not target_seq:
                    st.warning("Please provide a target sequence.")
                else:
                    aligner = Align.PairwiseAligner()
                    alignments = aligner.align(clean_seq, target_seq)
                    best_alignment = alignments[0]
                    
                    st.metric("Alignment Score", best_alignment.score)
                    with st.expander("View Alignment Data", expanded=True):
                        st.text(str(best_alignment))

        # --- TAB 5: Primer Analysis ---
        with tab_primer:
            st.markdown("### Primer Parameters")
            st.write("Calculate melting temperature (Tm) and GC content for custom primers.")
            
            if seq_type == "DNA":
                primer_input = st.text_input("Enter Primer Sequence (5' to 3')").upper().strip()
                if primer_input:
                    if set(primer_input).issubset(set("ATCGN")):
                        try:
                            tm_nn = mt.Tm_NN(Seq(primer_input))
                            primer_gc = gc_fraction(Seq(primer_input)) * 100
                            
                            p_col1, p_col2, p_col3 = st.columns(3)
                            p_col1.metric("Primer Length", f"{len(primer_input)} bp")
                            p_col2.metric("GC Content", f"{primer_gc:.1f}%")
                            p_col3.metric("Melting Temp (NN)", f"{tm_nn:.1f} °C")
                            
                            if primer_input in clean_seq:
                                st.success("Exact match found in primary sequence.")
                            elif str(Seq(primer_input).reverse_complement()) in clean_seq:
                                st.success("Exact match found in reverse complement of primary sequence.")
                            else:
                                st.info("No exact match found in primary sequence.")
                                
                        except Exception as e:
                            st.error(f"Calculation error: {e}")
                    else:
                        st.error("Primer contains invalid characters.")
            else:
                st.warning("Primer analysis requires a DNA sequence as the primary input.")
