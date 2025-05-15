import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Define a common font size
font_size_labels = 25
font_size_ticks = 22
font_size_legend = 23


import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# Load the custom Palatino font
font_path = "palatino.ttf"  # Update this with the correct path
palatino_font = fm.FontProperties(fname=font_path)

# Set the custom font as the default
plt.rcParams['font.family'] = palatino_font.get_name()


palatino_font_ticks = palatino_font.copy()
palatino_font_ticks.set_size(font_size_ticks)

palatino_font_labels = palatino_font.copy()
palatino_font_labels.set_size(font_size_labels)

palatino_font_legend = palatino_font.copy()
palatino_font_legend.set_size(font_size_legend)

# Read the CSV files
gc_data = pd.read_csv('plotGC.csv')
qa_data = pd.read_csv('plotQA.csv')
qk_data = pd.read_csv('plotQK.csv')
mc_data = pd.read_csv('plotMC.csv')

# Convert QPU access time from microseconds to milliseconds
gc_data['qpuaccesstime'] /= 1000
qa_data['qpuaccesstime'] /= 1000
qk_data['qpuaccesstime'] /= 1000
mc_data['qpuaccesstime'] /= 1000

# Calculate averages and variances for each numvar
gc_stats = gc_data.groupby('numvar')['qpuaccesstime'].agg(['mean', 'std']).reset_index()
qa_stats = qa_data.groupby('numvar')['qpuaccesstime'].agg(['mean', 'std']).reset_index()
qk_stats = qk_data.groupby('numvar')['qpuaccesstime'].agg(['mean', 'std']).reset_index()
mc_stats = mc_data.groupby('numvar')['qpuaccesstime'].agg(['mean', 'std']).reset_index()

# Replace any NaN in std with 0 (for groups with a single value)
gc_stats['std'] = gc_stats['std'].fillna(0)
qa_stats['std'] = qa_stats['std'].fillna(0)
qk_stats['std'] = qk_stats['std'].fillna(0)
mc_stats['std'] = mc_stats['std'].fillna(0)

# Create figure
plt.figure(figsize=(14, 14))

# Get the current axis
ax = plt.gca()

# Plot averages with lines only (no markers)
plt.plot(gc_stats['numvar'], gc_stats['mean'], label='GCP', color='#ca0020', linewidth=2)
plt.plot(qa_stats['numvar'], qa_stats['mean'], label='QAP', color='#f4a582', linewidth=2)
plt.plot(qk_stats['numvar'], qk_stats['mean'], label='QKP', color='#92c5de', linewidth=2)
plt.plot(mc_stats['numvar'], mc_stats['mean'], label='MCP', color='#0571b0', linewidth=2)

# Customize the plot
plt.xlabel('Number of Variables', fontproperties=palatino_font_labels)
plt.ylabel('QPU Access Time (ms)', fontproperties=palatino_font_labels)

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

ax.set_box_aspect(1)  # Keeps the plot area square
plt.subplots_adjust(left=0.15, right=0.85, top=0.85, bottom=0.3)

# Adjust layout to prevent label cutoff
#plt.tight_layout()

# Save the plot
plt.savefig('QPUAccessTime_MeanOnly.pdf', format='pdf', bbox_inches='tight')

# Show the plot
plt.show()
