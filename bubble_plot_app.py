import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import io

def round_to_nearest_ten(values):
    return [int(round(v / 10.0) * 10) for v in values]

def generate_bubble_plot(df, pathway_col, top_n, bubble_scale, cmap_choice, show_grid, plot_title):
    # Convert % to decimal for gene ratio
    df["GeneRatio"] = df["%"] / 100

    # Select top N based on p-value
    df = df.nsmallest(top_n, "PValue")
    df = df.sort_values(by="GeneRatio", ascending=True)

    # Color by -log10(p-value)
    df["neg_log_pval"] = -np.log10(df["PValue"])

    # Create the bubble plot
    fig, ax = plt.subplots(figsize=(12, 18))
    sns.set_style("whitegrid" if show_grid else "white")

    scatter = ax.scatter(
        df["GeneRatio"],
        df[pathway_col],
        s=df["Count"] * bubble_scale,
        c=df["neg_log_pval"],
        cmap=cmap_choice,
        edgecolors="black",
        alpha=0.8
    )

    if show_grid:
        ax.grid(True, linestyle="--", alpha=0.7)

    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.3, aspect=20)
    cbar.set_label("-log10 (p-value)", fontsize=14, fontweight="bold")

    # Gene Count legend
    min_count, max_count = df["Count"].min(), df["Count"].max()
    legend_sizes = np.linspace(min_count, max_count, num=3) if min_count != max_count else [min_count]
    legend_sizes = round_to_nearest_ten(legend_sizes)
    legend_sizes = sorted(set(legend_sizes))

    legend_handles = [ax.scatter([], [], s=size * bubble_scale * 0.8, color="gray", alpha=0.6) for size in legend_sizes]
    legend_labels = [f"{size}" for size in legend_sizes]

    ax.legend(handles=legend_handles, 
              labels=legend_labels, 
              title="Gene Counts",
              title_fontsize=12, 
              loc="upper left", 
              bbox_to_anchor=(1.05, 1), 
              fontsize=10, 
              frameon=True, 
              markerscale=1.0, 
              borderpad=1.5, 
              labelspacing=1.5, 
              handletextpad=1, 
              borderaxespad=1)

    # Labels and title
    ax.set_xlabel("Gene Ratio", fontsize=14, fontweight='bold')
    ax.set_ylabel("", fontsize=14, fontweight='bold')
    ax.set_title(plot_title, fontsize=14, fontweight='bold')
    ax.tick_params(axis='x', labelsize=14, width=2)
    ax.tick_params(axis='y', labelsize=16, width=2)

    return fig

# Streamlit UI
st.title("KEGG Pathway Bubble Plot Generator")
st.write("Upload a KEGG analysis Excel file to generate a bubble plot.")

uploaded_file = st.file_uploader("Upload Excel file", type=["xlsx"])

if uploaded_file:
    df = pd.read_excel(uploaded_file)
    st.write("Preview of uploaded data:", df.head())

    # Pathway column selector
    pathway_col = st.selectbox("Select the pathway column", df.columns, index=df.columns.get_loc("KEGGpathways") if "KEGGpathways" in df.columns else 0)

    top_n = st.slider("Select number of top pathways", 10, 50, 30)
    bubble_scale = st.slider("Adjust bubble scale", 10, 100, 20)
    cmap_choice = st.selectbox("Choose color map", ["RdBu_r", "viridis", "plasma", "coolwarm", "magma"])
    show_grid = st.checkbox("Show grid lines", True)
    plot_title = st.text_input("Enter plot title", "KEGG Pathway Analysis")

    st.write(f"**Current Figure Title:** {plot_title}")

    # Generate plot
    fig = generate_bubble_plot(df, pathway_col, top_n, bubble_scale, cmap_choice, show_grid, plot_title)
    st.pyplot(fig)

    # Download PDF
    buf = io.BytesIO()
    fig.savefig(buf, format="pdf", bbox_inches='tight')
    buf.seek(0)
    st.download_button("Download Plot as PDF", buf, file_name="KEGG_BubblePlot.pdf", mime="application/pdf")
