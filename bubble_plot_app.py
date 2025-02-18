import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import io

def generate_bubble_plot(df, top_n, bubble_scale, cmap_choice, show_grid):
    df["neg_log_pval"] = -np.log10(df["PValue"])
    df = df.nsmallest(top_n, "PValue")
    df = df.sort_values(by="neg_log_pval", ascending=True)

    plt.figure(figsize=(12, 18))
    plt.gca().yaxis.set_tick_params(pad=2)
    sns.set_style("whitegrid" if show_grid else "white")

    scatter = plt.scatter(
        df["neg_log_pval"], 
        df["KEGGpathways"], 
        s=df["Count"] * bubble_scale, 
        c=df["neg_log_pval"], 
        cmap=cmap_choice, 
        edgecolors="black", 
        alpha=0.8
    )

    if show_grid:
        plt.grid(True, linestyle="--", alpha=0.7)
    
    cbar = plt.colorbar(scatter, shrink=0.3, aspect=20)
    cbar.set_label("-log10 (p-value)", fontsize=14, fontweight="bold")

    legend_sizes = [10, 30, 50]
    for size in legend_sizes:
        plt.scatter([], [], s=size * bubble_scale, color="gray", alpha=0.6, label=f"Gene Count: {size}")
    
    plt.legend(title="Gene Counts", loc="upper right", bbox_to_anchor=(1.25, 1), fontsize=10)
    plt.xlabel("-log10 (p-value)", fontsize=14, fontweight='bold')
    plt.ylabel("", fontsize=14, fontweight='bold')
    plt.title("KEGG Upregulated Genes", fontsize=14, fontweight='bold')
    plt.xticks(fontsize=14, fontweight='bold')
    plt.yticks(fontsize=16, fontweight='bold')

    buf = io.BytesIO()
    plt.savefig(buf, format="pdf", bbox_inches='tight')
    plt.close()
    buf.seek(0)
    return buf

st.title("KEGG Bubble Plot Generator")
st.write("Upload an Excel file to generate a KEGG bubble plot.")

uploaded_file = st.file_uploader("Upload Excel file", type=["xlsx"])

if uploaded_file:
    df = pd.read_excel(uploaded_file)
    st.write("Preview of uploaded data:", df.head())
    
    top_n = st.slider("Select number of top pathways", 10, 50, 30)
    bubble_scale = st.slider("Adjust bubble scale", 10, 100, 20)
    cmap_choice = st.selectbox("Choose color map", ["RdBu_r", "viridis", "plasma", "coolwarm", "magma"])
    show_grid = st.checkbox("Show grid lines", True)
    
    if st.button("Generate Plot"):
        plot_buf = generate_bubble_plot(df, top_n, bubble_scale, cmap_choice, show_grid)
        st.download_button("Download Plot as PDF", plot_buf, file_name="KEGG_BubblePlot.pdf", mime="application/pdf")
