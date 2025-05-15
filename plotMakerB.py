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

# Find all unique numvar values
all_vars = sorted(set(gc_stats['numvar']).union(qa_stats['numvar'])
                                       .union(qk_stats['numvar'])
                                       .union(mc_stats['numvar']))

# Create figure with better font size and dimensions
plt.figure(figsize=(14, 10))

# Define a common font size
font_size_title = 18
font_size_labels = 18
font_size_ticks = 18
font_size_legend = 20

# Set the width of each bar and positions
bar_width = 0.2
x = np.arange(len(all_vars))

# Function to get values for specific variables
def get_values_for_vars(stats_df, all_vars):
    values = []
    for var in all_vars:
        value = stats_df[stats_df['numvar'] == var]['numqubit'].values
        values.append(value[0] if len(value) > 0 else 0)
    return values

# Get values for each dataset
gc_values = get_values_for_vars(gc_stats, all_vars)
qa_values = get_values_for_vars(qa_stats, all_vars)
qk_values = get_values_for_vars(qk_stats, all_vars)
mc_values = get_values_for_vars(mc_stats, all_vars)

# Create bars
plt.bar(x - bar_width*1.5, gc_values, width=bar_width, label='Graph Coloring', color='#ca0020')
plt.bar(x - bar_width*0.5, qa_values, width=bar_width, label='Quadratic Assignment', color='#f4a582')
plt.bar(x + bar_width*0.5, qk_values, width=bar_width, label='Quadratic Knapsack', color='#92c5de')
plt.bar(x + bar_width*1.5, mc_values, width=bar_width, label='Maximum Cut', color='#0571b0')

# Customize the plot
plt.xlabel('Number of Variables', fontsize=font_size_labels)
plt.ylabel('Number of Qubits', fontsize=font_size_labels)
plt.xticks(x, all_vars, fontsize=font_size_ticks)
plt.yticks(fontsize=font_size_ticks)
plt.grid(True, linestyle='--', alpha=0.7)

# Adjust legend to prevent overlap
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=font_size_legend)

# Adjust layout to prevent label cutoff
plt.tight_layout()

plt.savefig('NumberOfQubitsBar.pdf', format='pdf', bbox_inches='tight')

# Show the plot
plt.show()