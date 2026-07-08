import streamlit as st
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, molecular_weight
import pandas as pd

# 1. Page Configuration (Must be FIRST)
st.set_page_config(page_title="Zen Sequence Analyzer", page_icon="🌿", layout="wide")

# --- SIDEBAR: Keep inputs hidden away from the main canvas ---
with st.sidebar:
    st.title("🌿 Settings & Input")
    st.info("Paste your sequence below. The dashboard will update automatically.")
    
    # Input area is now tucked away
    dna_input = st.text_area("DNA Sequence (5' to 3')", "ATGCGATCGATCGATCGATCGATC", height=200)
    clean_dna = dna_input.strip().replace("\n", "").replace(" ", "").upper()
    
    st.markdown("---")
    st.caption("Built with Biopython & Streamlit")

# --- MAIN PAGE: Dedicated entirely to clean results ---
st.title("Sequence Analysis Dashboard")

if not clean_dna:
    st.warning("👈 Please enter a sequence in the sidebar to begin.")
else:
    valid_nucleotides = set("ATCGN")
    
    # Validation Check
    if not set(clean_dna).issubset(valid_nucleotides):
        st.error("⚠️ Invalid characters detected. Please input only A, T, C, G, or N.")
    else:
        # Core processing
        my_seq = Seq(clean_dna)
        
        # 2. TABS: Chunking information to prevent overwhelming the user
        tab_overview, tab_visuals, tab_biology = st.tabs([
            "📊 Quick Overview", 
            "📈 Visualizations", 
            "🧬 Translation & Downstream"
        ])
        
        # --- TAB 1: The Quick Numbers ---
        with tab_overview:
            st.markdown("### Core Metrics")
            # Creating 3 clean columns for KPI metrics
            col1, col2, col3 = st.columns(3)
            
            with col1:
                st.metric(label="Total Length", value=f"{len(clean_dna):,} bp")
            with col2:
                gc_content = gc_fraction(my_seq) * 100
                st.metric(label="GC Content", value=f"{gc_content:.1f}%")
            with col3:
                at_content = 100 - gc_content
                st.metric(label="AT Content", value=f"{at_content:.1f}%")
            
            st.markdown("---")
            # 3. EXPANDERS: Hide the ugly raw data to keep the vibe comfortable
            with st.expander("👀 View Raw Cleaned Sequence"):
                st.code(clean_dna, language="text")

        # --- TAB 2: The Charts ---
        with tab_visuals:
            st.markdown("### Nucleotide Distribution")
            
            # Calculate and display chart
            counts = {"A": clean_dna.count("A"), "T": clean_dna.count("T"), 
                      "C": clean_dna.count("C"), "G": clean_dna.count("G")}
            chart_data = pd.DataFrame(list(counts.values()), index=list(counts.keys()), columns=["Count"])
            
            # A nice, wide bar chart that breathes
            st.bar_chart(chart_data, height=300)

        # --- TAB 3: Advanced Biology ---
        with tab_biology:
            mrna = my_seq.transcribe()
            protein = my_seq.translate()
            
            st.markdown("### Biological Conversions")
            
            # Using columns inside the tab for balance
            bio_col1, bio_col2 = st.columns(2)
            
            with bio_col1:
                st.markdown("**Transcription (mRNA)**")
                with st.expander("View mRNA Sequence", expanded=True):
                    st.code(str(mrna), language="text")
            
            with bio_col2:
                st.markdown("**Translation (Protein)**")
                with st.expander("View Amino Acid Sequence", expanded=True):
                    st.code(str(protein), language="text")
                
                # Molecular weight tucked neatly at the bottom
                try:
                    clean_protein = str(protein).replace("*", "")
                    if clean_protein:
                        prot_weight = molecular_weight(clean_protein, seq_type="protein")
                        st.success(f"⚖️ **Estimated Mass:** {prot_weight:,.2f} Da")
                except Exception:
                    st.warning("Mass calculation unavailable for this sequence.")
