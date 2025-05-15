import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# Load the custom Palatino font
font_path = "palatino.ttf"  # Update this with the correct path
palatino_font = fm.FontProperties(fname=font_path)

# Set the custom font as the default
plt.rcParams['font.family'] = palatino_font.get_name()

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Read the CSV files
gc_data = pd.read_csv('plotGC.csv') 
qa_data = pd.read_csv('plotQA.csv')
qk_data = pd.read_csv('plotQK.csv')
mc_data = pd.read_csv('plotMC.csv')

# Calculate averages for each numvar
gc_stats = gc_data.groupby('numvar')['numqubit'].mean().reset_index()
qa_stats = qa_data.groupby('numvar')['numqubit'].mean().reset_index()
qk_stats = qk_data.groupby('numvar')['numqubit'].mean().reset_index()
mc_stats = mc_data.groupby('numvar')['numqubit'].mean().reset_index()

# Create figure with original dimensions
plt.figure(figsize=(14, 14))

# Get the current axis
ax = plt.gca()

# Define font sizes
font_size_labels = 25
font_size_ticks = 22
font_size_legend = 23

palatino_font_ticks = palatino_font.copy()
palatino_font_ticks.set_size(font_size_ticks)

palatino_font_labels = palatino_font.copy()
palatino_font_labels.set_size(font_size_labels)

palatino_font_legend = palatino_font.copy()
palatino_font_legend.set_size(font_size_legend)

# Plot averages with lines only (no markers)
plt.plot(gc_stats['numvar'], gc_stats['numqubit'], label='GCP', color='#ca0020', linewidth=2)
plt.plot(qa_stats['numvar'], qa_stats['numqubit'], label='QAP', color='#f4a582', linewidth=2)
plt.plot(qk_stats['numvar'], qk_stats['numqubit'], label='QKP', color='#92c5de', linewidth=2)
plt.plot(mc_stats['numvar'], mc_stats['numqubit'], label='MCP', color='#0571b0', linewidth=2)

plt.xlabel('Number of Variables', fontproperties=palatino_font_labels)
plt.ylabel('Number of Qubits', fontproperties=palatino_font_labels)

plt.xticks(fontproperties=palatino_font_ticks)
plt.yticks(fontproperties=palatino_font_ticks)

# Adjust legend with proper alignment
legend = plt.legend(
    loc='lower center', 
    bbox_to_anchor=(0.5, -0.35),  
    ncol=4,  
    frameon=True,  
    title="Optimization Problem",
    handletextpad=0.8,  # Space between line and text
    handleheight=2.0,  # Increase vertical alignment of lines
    handlelength=3.0,  # Make lines longer for better visibility
)

# Apply Palatino font to legend title and labels
legend.get_title().set_fontproperties(palatino_font_legend)
for text in legend.get_texts():
    text.set_fontproperties(palatino_font_legend)

plt.grid(True, linestyle='--', linewidth=0.7, alpha=0.6)

# Get the min/max values dynamically for both plots
y_min = min(gc_stats['numqubit'].min(), qa_stats['numqubit'].min(), 
            qk_stats['numqubit'].min(), mc_stats['numqubit'].min())

y_max = max(gc_stats['numqubit'].max(), qa_stats['numqubit'].max(), 
            qk_stats['numqubit'].max(), mc_stats['numqubit'].max())

# Apply to first plot
plt.ylim(y_min * 0.95, y_max * 1.05)  

ax.set_box_aspect(1)  # Keeps the plot area square
plt.subplots_adjust(left=0.15, right=0.85, top=0.85, bottom=0.3)

# Adjust layout to prevent label cutoff
#plt.tight_layout()

# Save the plot as a PDF
plt.savefig('NumberOfQubits_MeanOnly.pdf', format='pdf', bbox_inches='tight')

# Show the plot
plt.show()
