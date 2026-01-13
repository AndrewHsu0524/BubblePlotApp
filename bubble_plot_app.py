import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import io

def round_to_nearest_ten(values):
    return [int(round(v / 10.0) * 10) for v in values]

def parse_percent_series(s: pd.Series) -> pd.Series:
    if not pd.api.types.is_numeric_dtype(s):
        extracted = s.astype(str).str.extract(r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)")[0]
        s_num = pd.to_numeric(extracted, errors="coerce")
    else:
        s_num = pd.to_numeric(s, errors="coerce")

    valid = s_num.dropna()
    if len(valid) == 0:
        return s_num

    ratio = s_num.copy()
    if (valid > 1).mean() > 0.5:
        ratio = s_num / 100.0
    return ratio

def auto_pick_column(columns, candidates):
    for key in candidates:
        for c in columns:
            if key in c.lower():
                return c
    return None

def generate_bubble_plot(
    df,
    pathway_col,
    percent_col,
    sig_col,
    count_col,
    top_n,
    sig_cutoff,
    bubble_scale,
    cmap_choice,
    show_grid,
    plot_title,
    sort_mode,
    fig_w,
    fig_h,
    title_fs,
    y_fs,
    x_label_fs,
    x_tick_fs,
    transparent_bg,
):
    df = df.copy()

    df["GeneRatio"] = parse_percent_series(df[percent_col])
    df["Sig_num"] = pd.to_numeric(df[sig_col], errors="coerce")
    df["Count_num"] = pd.to_numeric(df[count_col], errors="coerce")

    df = df.dropna(subset=["GeneRatio", "Sig_num", "Count_num", pathway_col])
    if df.empty:
        raise ValueError("No valid rows after cleaning. Check selected columns / data format.")

    df = df[df["Sig_num"] > 0]
    if df.empty:
        raise ValueError("All significance values are <= 0. Cannot compute -log10.")

    if sig_cutoff is not None:
        df = df[df["Sig_num"] <= sig_cutoff]
        if df.empty:
            raise ValueError(f"No rows pass the significance cutoff (≤ {sig_cutoff}).")

    df = df.nsmallest(top_n, "Sig_num")

    if sort_mode == "GeneRatio":
        df = df.sort_values(by="GeneRatio", ascending=True)
    elif sort_mode == "Significance":
        df = df.sort_values(by="Sig_num", ascending=True)
    elif sort_mode == "Alphabetical":
        df = df.sort_values(by=pathway_col, ascending=True)

    df["neg_log_sig"] = -np.log10(df["Sig_num"])

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.set_style("whitegrid" if show_grid else "white")

    scatter = ax.scatter(
        df["GeneRatio"],
        df[pathway_col].astype(str),
        s=df["Count_num"] * bubble_scale,
        c=df["neg_log_sig"],
        cmap=cmap_choice,
        edgecolors="black",
        alpha=0.8,
    )

    if show_grid:
        ax.grid(True, linestyle="--", alpha=0.7)

    cbar = plt.colorbar(scatter, ax=ax, shrink=0.3, aspect=20)
    cbar.set_label("-log10 (Significance)", fontsize=max(8, int(x_label_fs)), fontweight="bold")

    min_count, max_count = df["Count_num"].min(), df["Count_num"].max()
    legend_sizes = np.linspace(min_count, max_count, num=3) if min_count != max_count else [min_count]
    legend_sizes = round_to_nearest_ten(legend_sizes)
    legend_sizes = sorted(set(int(x) for x in legend_sizes if not pd.isna(x)))
    if len(legend_sizes) == 0:
        legend_sizes = [10]

    legend_handles = [
        ax.scatter([], [], s=size * bubble_scale * 0.8, color="gray", alpha=0.6)
        for size in legend_sizes
    ]
    legend_labels = [f"{size}" for size in legend_sizes]

    ax.legend(
        handles=legend_handles,
        labels=legend_labels,
        title="Gene Counts",
        title_fontsize=max(8, int(y_fs)),
        loc="upper left",
        bbox_to_anchor=(1.05, 1),
        fontsize=max(8, int(y_fs) - 2),
        frameon=True,
    )

    ax.set_xlabel("Gene Ratio", fontsize=x_label_fs, fontweight="bold")
    ax.set_ylabel("", fontsize=y_fs, fontweight="bold")
    ax.set_title(plot_title, fontsize=title_fs, fontweight="bold")
    ax.tick_params(axis="x", labelsize=x_tick_fs, width=2)
    ax.tick_params(axis="y", labelsize=y_fs, width=2)

    if transparent_bg:
        fig.patch.set_alpha(0)
        ax.patch.set_alpha(0)

    return fig

# -------------------------
# App
# -------------------------
st.set_page_config(layout="wide")

st.title("KEGG Pathway Bubble Plot Generator")
st.write("Upload an Excel file. Adjust settings in the sidebar; preview stays centered on this page.")

uploaded_file = st.file_uploader("Upload Excel file", type=["xlsx"])

if uploaded_file:
    df = pd.read_excel(uploaded_file)

    # Full-width preview table (middle)
    st.subheader("Data preview")
    st.dataframe(df, use_container_width=True, height=260)
    st.caption(f"{len(df)} rows × {len(df.columns)} columns")
    st.markdown("---")

    # Auto-detect
    default_pathway = auto_pick_column(df.columns, ["kegg", "pathway", "term", "description"])
    default_percent = auto_pick_column(df.columns, ["%", "percent", "percentage", "gene ratio", "generatio", "ratio"])
    default_sig = auto_pick_column(df.columns, ["fdr", "qvalue", "q-value", "adj", "adjust", "padj", "p.adjust", "pvalue", "p-value"])
    default_count = auto_pick_column(df.columns, ["count", "genes", "gene count", "hits"])

    # -------------------------
    # Sidebar controls (scrollable)
    # -------------------------
    st.sidebar.header("Settings")

    with st.sidebar.expander("1) Columns", expanded=True):
        pathway_col = st.selectbox(
            "Pathway / Term column",
            df.columns,
            index=df.columns.get_loc(default_pathway) if default_pathway in df.columns else 0,
        )
        percent_col = st.selectbox(
            "Percent / Gene ratio column",
            df.columns,
            index=df.columns.get_loc(default_percent) if default_percent in df.columns else 0,
        )
        sig_col = st.selectbox(
            "Significance (FDR / adjusted p-value / p-value)",
            df.columns,
            index=df.columns.get_loc(default_sig) if default_sig in df.columns else 0,
            help="Smaller = more significant.",
        )
        count_col = st.selectbox(
            "Gene Count column",
            df.columns,
            index=df.columns.get_loc(default_count) if default_count in df.columns else 0,
        )

    with st.sidebar.expander("2) Filters & Sorting", expanded=True):
        top_n = st.slider("Top pathways (by Significance)", 10, 50, 30)
        enable_cutoff = st.checkbox("Enable significance cutoff", True)
        sig_cutoff = None
        if enable_cutoff:
            sig_cutoff = st.number_input("Significance cutoff (≤)", value=0.05, min_value=0.0, format="%.6f")
        sort_mode = st.selectbox("Sort pathways by", ["GeneRatio", "Significance", "Alphabetical"], index=0)

    with st.sidebar.expander("3) Style", expanded=False):
        plot_title = st.text_input("Plot title", "KEGG Pathway Analysis")
        show_grid = st.checkbox("Show grid lines", True)
        cmap_choice = st.selectbox("Color map", ["RdBu_r", "viridis", "plasma", "coolwarm", "magma"])
        bubble_scale = st.slider("Bubble scale", 10, 100, 20)

        st.markdown("**Figure size**")
        fig_w = st.slider("Width", 6, 24, 12)
        fig_h = st.slider("Height", 6, 36, 18)

        st.markdown("**Font sizes**")
        title_fs = st.slider("Title", 10, 30, 16)
        y_fs = st.slider("Y-axis (pathways)", 8, 30, 16)
        x_label_fs = st.slider("X-axis label", 8, 30, 14)
        x_tick_fs = st.slider("X-axis ticks", 6, 24, 12)

    with st.sidebar.expander("4) Export", expanded=False):
        export_format = st.selectbox("Export format", ["pdf", "png", "svg"], index=0)
        transparent_bg = st.checkbox("Transparent background", False)
        png_dpi = st.number_input("PNG DPI", value=300, min_value=72, max_value=1200, step=50)

    # Provide safe defaults if user doesn't open style/export
    # (Streamlit requires variables exist before use)
    if "plot_title" not in locals():
        plot_title = "KEGG Pathway Analysis"
        show_grid = True
        cmap_choice = "RdBu_r"
        bubble_scale = 20
        fig_w, fig_h = 12, 18
        title_fs, y_fs, x_label_fs, x_tick_fs = 16, 16, 14, 12
    if "export_format" not in locals():
        export_format = "pdf"
        transparent_bg = False
        png_dpi = 300

    # -------------------------
    # Main area: plot stays centered
    # -------------------------
    st.subheader("Plot preview")

    try:
        fig = generate_bubble_plot(
            df=df,
            pathway_col=pathway_col,
            percent_col=percent_col,
            sig_col=sig_col,
            count_col=count_col,
            top_n=top_n,
            sig_cutoff=sig_cutoff,
            bubble_scale=bubble_scale,
            cmap_choice=cmap_choice,
            show_grid=show_grid,
            plot_title=plot_title,
            sort_mode=sort_mode,
            fig_w=fig_w,
            fig_h=fig_h,
            title_fs=title_fs,
            y_fs=y_fs,
            x_label_fs=x_label_fs,
            x_tick_fs=x_tick_fs,
            transparent_bg=transparent_bg,
        )

        st.pyplot(fig, use_container_width=True)

        buf = io.BytesIO()
        save_kwargs = dict(bbox_inches="tight", transparent=transparent_bg)

        if export_format == "png":
            fig.savefig(buf, format="png", dpi=png_dpi, **save_kwargs)
            mime = "image/png"
            fname = "KEGG_BubblePlot.png"
        elif export_format == "svg":
            fig.savefig(buf, format="svg", **save_kwargs)
            mime = "image/svg+xml"
            fname = "KEGG_BubblePlot.svg"
        else:
            fig.savefig(buf, format="pdf", **save_kwargs)
            mime = "application/pdf"
            fname = "KEGG_BubblePlot.pdf"

        buf.seek(0)
        st.download_button(
            f"Download {export_format.upper()}",
            buf,
            file_name=fname,
            mime=mime,
            use_container_width=False,
        )

    except Exception as e:
        st.error(f"Could not generate plot: {e}")
        st.info("Tip: double-check the selected columns and ensure Significance > 0 and Count is numeric.")
