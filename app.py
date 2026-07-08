import streamlit as st
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, molecular_weight
import pandas as pd

# 1. Page Configuration (Must be FIRST)
# Removed emojis and set a clean title to match the wireframe
st.set_page_config(page_title="Biomolecule Sequence Analyser", layout="wide")

# Initialize session state for sequence history (from wireframe)
if "sequence_history" not in st.session_state:
    st.session_state.sequence_history = []

# --- SIDEBAR: Settings, Input, and History ---
with st.sidebar:
    st.title("Biomolecule Sequence Analyser")
    st.markdown("Welcome back. Ready to analyze a new sequence?")
    
    # 2. Sequence Type Selector (Added Protein and RNA)
    seq_type = st.selectbox("Sequence Type", ["DNA", "RNA", "Protein"])
    
    # Dynamic placeholder based on selected type
    default_seq = ""
    if seq_type == "DNA":
        default_seq = "ATGCGATCGATCGATCGATCGATC"
    elif seq_type == "RNA":
        default_seq = "AUGCGAUCGAUCGAUCGAUCGAUC"
    else:
        default_seq = "MKWVTFISLLFLFSSAYSRGVFRR" # Example protein
    
    # Input area
    raw_input = st.text_area("Paste your sequence here", default_seq, height=150)
    clean_seq = raw_input.strip().replace("\n", "").replace(" ", "").upper()
    
    # Analyze button to trigger history recording
    if st.button("Analyze Sequence"):
        if clean_seq and clean_seq not in st.session_state.sequence_history:
            st.session_state.sequence_history.insert(0, clean_seq) # Add newest to top
            st.session_state.sequence_history = st.session_state.sequence_history[:5] # Keep last 5
            
    st.markdown("---")
    
    # Display History
    st.subheader("Previously Analysed")
    if not st.session_state.sequence_history:
        st.caption("No history in this session.")
    else:
        for i, seq in enumerate(st.session_state.sequence_history):
            # Show a shortened preview of the sequence
            st.text(f"{i+1}. {seq[:15]}...")

# --- MAIN PAGE: Clean Results ---
if not clean_seq:
    st.info("Please enter a sequence in the sidebar to begin.")
else:
    # 3. Dynamic Validation Check based on type
    is_valid = True
    error_msg = ""
    
    if seq_type == "DNA":
        valid_chars = set("ATCGN")
        if not set(clean_seq).issubset(valid_chars):
            is_valid = False
            error_msg = "Invalid characters detected. Please input only A, T, C, G, or N."
    elif seq_type == "RNA":
        valid_chars = set("AUCGN")
        if not set(clean_seq).issubset(valid_chars):
            is_valid = False
            error_msg = "Invalid characters detected. Please input only A, U, C, G, or N."
    elif seq_type == "Protein":
        # Standard IUPAC amino acids plus asterisk for stop codons
        valid_chars = set("ACDEFGHIKLMNPQRSTVWY*")
        if not set(clean_seq).issubset(valid_chars):
            is_valid = False
            error_msg = "Invalid characters detected. Please input valid standard amino acids."

    if not is_valid:
        st.error(error_msg)
    else:
        # Core processing
        my_seq = Seq(clean_seq)
        
        # TABS: Chunking information
        tab_overview, tab_visuals, tab_biology = st.tabs([
            "Quick Overview", 
            "Visualizations", 
            "Downstream Analysis"
        ])
        
        # --- TAB 1: The Quick Numbers ---
        with tab_overview:
            st.markdown("### Core Metrics")
            
            # Metrics change depending on if it's a nucleic acid or protein
            if seq_type in ["DNA", "RNA"]:
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric(label="Total Length", value=f"{len(clean_seq):,} bp")
                with col2:
                    gc_content = gc_fraction(my_seq) * 100
                    st.metric(label="GC Content", value=f"{gc_content:.1f}%")
                with col3:
                    at_content = 100 - gc_content
                    st.metric(label="AT/AU Content", value=f"{at_content:.1f}%")
            else: # Protein
                col1, col2 = st.columns(2)
                with col1:
                    st.metric(label="Total Length", value=f"{len(clean_seq):,} AA")
                with col2:
                    try:
                        clean_protein = clean_seq.replace("*", "")
                        prot_weight = molecular_weight(clean_protein, seq_type="protein")
                        st.metric(label="Estimated Mass", value=f"{prot_weight:,.2f} Da")
                    except Exception:
                        st.metric(label="Estimated Mass", value="N/A")
            
            st.markdown("---")
            with st.expander("View Raw Cleaned Sequence"):
                st.code(clean_seq, language="text")

        # --- TAB 2: The Charts ---
        with tab_visuals:
            st.markdown(f"### {seq_type} Distribution")
            
            if seq_type == "DNA":
                counts = {"A": clean_seq.count("A"), "T": clean_seq.count("T"), 
                          "C": clean_seq.count("C"), "G": clean_seq.count("G")}
            elif seq_type == "RNA":
                counts = {"A": clean_seq.count("A"), "U": clean_seq.count("U"), 
                          "C": clean_seq.count("C"), "G": clean_seq.count("G")}
            else: # Protein
                counts = {}
                for aa in set(clean_seq):
                    if aa in valid_chars and aa != "*":
                        counts[aa] = clean_seq.count(aa)
                counts = dict(sorted(counts.items())) # Sort alphabetically
                
            chart_data = pd.DataFrame(list(counts.values()), index=list(counts.keys()), columns=["Count"])
            st.bar_chart(chart_data, height=300)

        # --- TAB 3: Advanced Biology ---
        with tab_biology:
            if seq_type == "DNA":
                st.markdown("### Biological Conversions")
                mrna = my_seq.transcribe()
                protein = my_seq.translate()
                
                bio_col1, bio_col2 = st.columns(2)
                
                with bio_col1:
                    st.markdown("**Transcription (mRNA)**")
                    with st.expander("View mRNA Sequence", expanded=True):
                        st.code(str(mrna), language="text")
                
                with bio_col2:
                    st.markdown("**Translation (Protein)**")
                    with st.expander("View Amino Acid Sequence", expanded=True):
                        st.code(str(protein), language="text")
                    
                    try:
                        clean_protein = str(protein).replace("*", "")
                        if clean_protein:
                            prot_weight = molecular_weight(clean_protein, seq_type="protein")
                            st.info(f"Estimated Protein Mass: {prot_weight:,.2f} Da")
                    except Exception:
                        st.warning("Mass calculation unavailable for this sequence.")
            
            elif seq_type == "RNA":
                st.markdown("### Biological Conversions")
                protein = my_seq.translate()
                st.markdown("**Translation (Protein)**")
                with st.expander("View Amino Acid Sequence", expanded=True):
                    st.code(str(protein), language="text")
                try:
                    clean_protein = str(protein).replace("*", "")
                    if clean_protein:
                        prot_weight = molecular_weight(clean_protein, seq_type="protein")
                        st.info(f"Estimated Protein Mass: {prot_weight:,.2f} Da")
                except Exception:
                    st.warning("Mass calculation unavailable for this sequence.")
            
            elif seq_type == "Protein":
                st.markdown("### Protein Properties")
                st.write("This sequence is already a protein. Reverse translation to a specific nucleotide sequence cannot be accurately predicted due to codon degeneracy.")
                try:
                    clean_protein = clean_seq.replace("*", "")
                    if clean_protein:
                        prot_weight = molecular_weight(clean_protein, seq_type="protein")
                        st.info(f"Estimated Mass: {prot_weight:,.2f} Daltons")
                except Exception:
                    st.warning("Mass calculation unavailable for this sequence.")
