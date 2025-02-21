**Pathway Analysis Bubble Plot Generator**

This is a Streamlit-based web application that allows users to generate interactive bubble plots for pathway analysis. 
Users can upload an Excel file containing pathway analysis results and dynamically adjust plot parameters with live previews.

**Features**
  - Upload an Excel file containing pathway data
  - Select the pathway column dynamically
  - Adjust visualization settings (number of top pathways, bubble scale, color map, grid visibility, etc.)
  - Live preview of the plot as users make changes
  - Download the final plot as a PDF

**Running the App on Streamlit Cloud**
  You can run this app directly on Streamlit Cloud without any local installation. Follow these steps:
  1. Fork this repository or upload the bubble_plot_app.py file to your GitHub repository.
  2. Go to Streamlit Community Cloud and sign in with GitHub.
  3. Click "New App" and select your repository.
  4. Choose bubble_plot_app.py as the main script.
  5. Click "Deploy" and Streamlit will host your app online!
  6. Share the Streamlit URL with others to allow them to use the tool without installation.

**Local Installation**
  1. Clone the repository:
     git clone https://github.com/yourusername/PathwayBubblePlot.git
     cd PathwayBubblePlot
  2. Install dependencies:
     pip install -r requirements.txt
  3. Run the application:
     streamlit run bubble_plot_app.py

**Usage**
  1. Upload an Excel file with pathway data (including columns for pathway names, p-values, and gene counts).
  2. Customize the visualization with sliders and dropdowns.
  3. Preview the plot live as you make adjustments.
  4. Download the final figure as a high-quality PDF.

**Example Input Data Format**
  The uploaded Excel file should contain the following columns:
  | KEGGpathways | PValue  | Count |
  |-------------|---------|--------|
  | Pathway A   | 0.0023  | 25     |
  | Pathway B   | 0.0081  | 40     |
  | Pathway C   | 0.0150  | 60     |

**Dependencies**
  - Python 3.8+
  - Streamlit
  - Pandas
  - Matplotlib
  - Seaborn
  - NumPy
  - Openpyxl
