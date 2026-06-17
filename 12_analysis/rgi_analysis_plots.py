import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import glob
import os
import warnings
warnings.filterwarnings('ignore')

print("Libraries loaded successfully!")

# ── Load and Combine All TSV Files ───────────────────────────
files_list = glob.glob('*_combined_rgi_bins.tsv')

dfs = []
for f in files_list:
    sample = os.path.basename(f).replace('_combined_rgi_bins.tsv', '')
    df = pd.read_csv(f, sep='\t')
    df['Sample'] = sample
    dfs.append(df)

all_df = pd.concat(dfs, ignore_index=True)
print(f"Loaded {len(files_list)} files — {all_df.shape[0]} total rows")

# ── Filter to Strict + Perfect hits only ─────────────────────
filt = all_df[all_df['Cut_Off'].isin(['Strict', 'Perfect'])].copy()
print(f"After filtering (Strict + Perfect): {filt.shape[0]} hits")

# ── Extract Genus from Taxonomy column ───────────────────────
def get_genus(tax):
    try:
        for p in str(tax).split(';'):
            if p.strip().startswith('g__'):
                return p.strip().replace('g__', '')
        return 'Unknown'
    except:
        return 'Unknown'

filt['Genus'] = filt['Taxonomy'].apply(get_genus)

# ── Sort samples numerically ──────────────────────────────────
def sort_key(s):
    try:
        return int(s)
    except:
        return 9999

sample_order = sorted(filt['Sample'].unique(), key=sort_key)
print(f"Samples found: {sample_order}")


# ════════════════════════════════════════════════════════════
#  PLOT 1 — AMR Gene Hit Count per Sample (Bar Chart)
# ════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(12, 5))

counts = filt.groupby('Sample').size().reindex(sample_order)
colors = plt.cm.viridis(np.linspace(0.2, 0.85, len(counts)))

bars = ax.bar(sample_order, counts.values, color=colors,
              edgecolor='white', linewidth=0.8, width=0.65)

for bar, v in zip(bars, counts.values):
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1.5,
            str(v), ha='center', va='bottom', fontsize=9, fontweight='bold')

ax.set_xlabel('Sample', fontsize=12)
ax.set_ylabel('Number of AMR Hits', fontsize=12)
ax.set_title('AMR Gene Hits per Sample (Strict + Perfect)', fontsize=14,
             fontweight='bold', pad=12)
ax.set_facecolor('#f8f8f8')
fig.patch.set_facecolor('white')
ax.spines[['top', 'right']].set_visible(False)
ax.set_ylim(0, counts.max() * 1.12)

plt.tight_layout()
plt.savefig('01_amr_counts_per_sample.png', dpi=180, bbox_inches='tight')
plt.show()
print("Plot 1 saved: 01_amr_counts_per_sample.png")


# ════════════════════════════════════════════════════════════
#  PLOT 2 — Drug Class Distribution per Sample (Stacked Bar)
# ════════════════════════════════════════════════════════════

dc_sample = (filt.groupby(['Sample', 'Drug_Class'])
             .size().unstack(fill_value=0).reindex(sample_order))
top_dc = dc_sample.sum().nlargest(12).index
dc_plot = dc_sample[top_dc]

cmap_dc = plt.cm.get_cmap('tab20', len(top_dc))
colors_dc = [cmap_dc(i) for i in range(len(top_dc))]

fig, ax = plt.subplots(figsize=(16, 7))
bottom = np.zeros(len(sample_order))

for i, dc in enumerate(top_dc):
    vals = dc_plot[dc].values
    ax.bar(sample_order, vals, bottom=bottom, label=dc,
           color=colors_dc[i], edgecolor='white', linewidth=0.6, width=0.7)
    bottom += vals

ax.set_xlabel('Sample', fontsize=12)
ax.set_ylabel('AMR Gene Count', fontsize=12)
ax.set_title('Drug Class Distribution per Sample', fontsize=14,
             fontweight='bold', pad=12)
ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=8.5,
          framealpha=0.9, title='Drug Class', title_fontsize=9)
ax.set_facecolor('#f8f8f8')
fig.patch.set_facecolor('white')
ax.spines[['top', 'right']].set_visible(False)

plt.tight_layout()
plt.savefig('02_drug_class_stacked_bar.png', dpi=180, bbox_inches='tight')
plt.show()
print("Plot 2 saved: 02_drug_class_stacked_bar.png")


# ════════════════════════════════════════════════════════════
#  PLOT 3 — Resistance Mechanism per Sample (Stacked Bar)
# ════════════════════════════════════════════════════════════

mech_sample = (filt.groupby(['Sample', 'Resistance_Mechanism'])
               .size().unstack(fill_value=0).reindex(sample_order))
mechs = mech_sample.columns.tolist()
cmap_m = plt.cm.get_cmap('Set2', len(mechs))

fig, ax = plt.subplots(figsize=(16, 7))
bottom = np.zeros(len(sample_order))

for i, m in enumerate(mechs):
    vals = mech_sample[m].values
    ax.bar(sample_order, vals, bottom=bottom, label=m,
           color=cmap_m(i), edgecolor='white', linewidth=0.6, width=0.7)
    bottom += vals

ax.set_xlabel('Sample', fontsize=12)
ax.set_ylabel('AMR Gene Count', fontsize=12)
ax.set_title('Resistance Mechanism Distribution per Sample', fontsize=14,
             fontweight='bold', pad=12)
ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=8.5,
          framealpha=0.9, title='Mechanism', title_fontsize=9)
ax.set_facecolor('#f8f8f8')
fig.patch.set_facecolor('white')
ax.spines[['top', 'right']].set_visible(False)

plt.tight_layout()
plt.savefig('03_resistance_mechanism_stacked.png', dpi=180, bbox_inches='tight')
plt.show()
print("Plot 3 saved: 03_resistance_mechanism_stacked.png")


# ════════════════════════════════════════════════════════════
#  PLOT 4 — Overall Resistance Mechanism (Pie Chart)
# ════════════════════════════════════════════════════════════

mech_counts = filt['Resistance_Mechanism'].value_counts()
colors_pie = plt.cm.Set3(np.linspace(0, 1, len(mech_counts)))

fig, ax = plt.subplots(figsize=(9, 7))
wedges, texts, autotexts = ax.pie(
    mech_counts.values,
    labels=None,
    autopct='%1.1f%%',
    colors=colors_pie,
    startangle=140,
    pctdistance=0.75,
    wedgeprops=dict(edgecolor='white', linewidth=1.5)
)
for at in autotexts:
    at.set_fontsize(9)

ax.legend(wedges, mech_counts.index,
          loc='lower center', bbox_to_anchor=(0.5, -0.18),
          ncol=2, fontsize=9, framealpha=0.9,
          title='Resistance Mechanism')
ax.set_title('Overall Distribution of Resistance Mechanisms',
             fontsize=13, fontweight='bold', pad=15)
fig.patch.set_facecolor('white')

plt.tight_layout()
plt.savefig('04_resistance_mechanism_pie.png', dpi=180, bbox_inches='tight')
plt.show()
print("Plot 4 saved: 04_resistance_mechanism_pie.png")


# ════════════════════════════════════════════════════════════
#  PLOT 5 — Heatmap: Drug Class x Sample
# ════════════════════════════════════════════════════════════

heat_dc = dc_plot.T
cmap_heat = LinearSegmentedColormap.from_list(
    'wbu', ['#ffffff', '#c6dbef', '#2171b5', '#08306b'])

fig, ax = plt.subplots(figsize=(13, 7))
sns.heatmap(heat_dc, ax=ax, cmap=cmap_heat,
            linewidths=0.5, linecolor='#dddddd',
            annot=True, fmt='d', annot_kws={'size': 8},
            cbar_kws={'label': 'AMR Hit Count', 'shrink': 0.7})

ax.set_xlabel('Sample', fontsize=11)
ax.set_ylabel('Drug Class', fontsize=11)
ax.set_title('Heatmap: Drug Class x Sample', fontsize=14,
             fontweight='bold', pad=12)
plt.xticks(rotation=45, ha='right', fontsize=9)
plt.yticks(fontsize=8.5)
fig.patch.set_facecolor('white')

plt.tight_layout()
plt.savefig('05_heatmap_drugclass_sample.png', dpi=180, bbox_inches='tight')
plt.show()
print("Plot 5 saved: 05_heatmap_drugclass_sample.png")


# ════════════════════════════════════════════════════════════
#  PLOT 6 — Heatmap: AMR Gene Family x Sample (Top 20)
# ════════════════════════════════════════════════════════════

fam_sample = (filt.groupby(['Sample', 'AMR_Gene_Family'])
              .size().unstack(fill_value=0).reindex(sample_order))
top_fam = fam_sample.sum().nlargest(20).index
fam_plot = fam_sample[top_fam].T

fam_plot.index = [
    i[:55] + '...' if len(i) > 55 else i for i in fam_plot.index
]

fig, ax = plt.subplots(figsize=(14, 9))
sns.heatmap(fam_plot, ax=ax, cmap='YlOrRd',
            linewidths=0.4, linecolor='#eeeeee',
            annot=True, fmt='d', annot_kws={'size': 7.5},
            cbar_kws={'label': 'AMR Hit Count', 'shrink': 0.65})

ax.set_xlabel('Sample', fontsize=11)
ax.set_ylabel('AMR Gene Family', fontsize=11)
ax.set_title('Heatmap: AMR Gene Family x Sample (Top 20)',
             fontsize=14, fontweight='bold', pad=12)
plt.xticks(rotation=45, ha='right', fontsize=9)
plt.yticks(fontsize=7.5)
fig.patch.set_facecolor('white')

plt.tight_layout()
plt.savefig('06_heatmap_genefamily_sample.png', dpi=180, bbox_inches='tight')
plt.show()
print("Plot 6 saved: 06_heatmap_genefamily_sample.png")


# ════════════════════════════════════════════════════════════
#  PLOT 7 — Top 20 AMR Gene Families (Horizontal Bar)
# ════════════════════════════════════════════════════════════

fam_counts = filt['AMR_Gene_Family'].value_counts().head(20)
short_labels = [l[:60] + '...' if len(l) > 60 else l for l in fam_counts.index]
colors_bar = plt.cm.plasma(np.linspace(0.15, 0.85, len(fam_counts)))

fig, ax = plt.subplots(figsize=(11, 8))
bars = ax.barh(short_labels[::-1], fam_counts.values[::-1],
               color=colors_bar[::-1], edgecolor='white', height=0.7)

for bar, v in zip(bars, fam_counts.values[::-1]):
    ax.text(v + 1, bar.get_y() + bar.get_height() / 2, str(v),
            va='center', fontsize=8.5, fontweight='bold')

ax.set_xlabel('Count', fontsize=11)
ax.set_title('Top 20 AMR Gene Families (All Samples)',
             fontsize=13, fontweight='bold', pad=12)
ax.set_facecolor('#f8f8f8')
fig.patch.set_facecolor('white')
ax.spines[['top', 'right']].set_visible(False)
ax.set_xlim(0, fam_counts.max() * 1.13)

plt.tight_layout()
plt.savefig('07_top_amr_gene_families.png', dpi=180, bbox_inches='tight')
plt.show()
print("Plot 7 saved: 07_top_amr_gene_families.png")


print("\nAll plots 1-7 generated successfully!")

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os
import warnings
warnings.filterwarnings('ignore')

print("Libraries loaded successfully!")

# ── Load and Combine All TSV Files ───────────────────────────
files_list = glob.glob('*_combined_rgi_bins.tsv')

dfs = []
for f in files_list:
    sample = os.path.basename(f).replace('_combined_rgi_bins.tsv', '')
    df = pd.read_csv(f, sep='\t')
    df['Sample'] = sample
    dfs.append(df)

all_df = pd.concat(dfs, ignore_index=True)

# ── Filter to Strict + Perfect hits only ─────────────────────
filt = all_df[all_df['Cut_Off'].isin(['Strict', 'Perfect'])].copy()
print(f"After filtering (Strict + Perfect): {filt.shape[0]} hits")

# ── Extract Genus from Taxonomy column ───────────────────────
def get_genus(tax):
    try:
        for p in str(tax).split(';'):
            if p.strip().startswith('g__'):
                return p.strip().replace('g__', '')
        return 'Unknown'
    except:
        return 'Unknown'

filt['Genus'] = filt['Taxonomy'].apply(get_genus)

# ── Sort samples numerically ──────────────────────────────────
def sort_key(s):
    try:
        return int(s)
    except:
        return 9999

sample_order = sorted(filt['Sample'].unique(), key=sort_key)


# ════════════════════════════════════════════════════════════
#  PLOT 8 — Top 15 Genera by AMR Gene Count (Bar Chart)
# ════════════════════════════════════════════════════════════

genus_counts = filt['Genus'].value_counts().head(15)
colors_g = plt.cm.tab20(np.linspace(0, 1, len(genus_counts)))

fig, ax = plt.subplots(figsize=(11, 6))
bars = ax.bar(genus_counts.index, genus_counts.values,
              color=colors_g, edgecolor='white', linewidth=0.8, width=0.7)

for bar, v in zip(bars, genus_counts.values):
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
            str(v), ha='center', va='bottom', fontsize=9, fontweight='bold')

ax.set_xlabel('Genus', fontsize=11)
ax.set_ylabel('AMR Gene Count', fontsize=11)
ax.set_title('Top 15 Genera by AMR Gene Count',
             fontsize=13, fontweight='bold', pad=12)
plt.xticks(rotation=40, ha='right', fontsize=9)
ax.set_facecolor('#f8f8f8')
fig.patch.set_facecolor('white')
ax.spines[['top', 'right']].set_visible(False)
ax.set_ylim(0, genus_counts.max() * 1.12)

plt.tight_layout()
plt.savefig('09_top_genera_amr.png', dpi=180, bbox_inches='tight')
plt.show()
print("Plot 9 saved: 09_top_genera_amr.png")


# ════════════════════════════════════════════════════════════
#  PLOT 9 — Genus x Resistance Mechanism Heatmap (Top 12)
# ════════════════════════════════════════════════════════════

top_genera = filt['Genus'].value_counts().head(12).index
gen_mech = (filt[filt['Genus'].isin(top_genera)]
            .groupby(['Genus', 'Resistance_Mechanism'])
            .size().unstack(fill_value=0))

fig, ax = plt.subplots(figsize=(12, 7))
sns.heatmap(gen_mech, ax=ax, cmap='Greens',
            linewidths=0.5, linecolor='#dddddd',
            annot=True, fmt='d', annot_kws={'size': 9},
            cbar_kws={'label': 'Count', 'shrink': 0.65})

ax.set_xlabel('Resistance Mechanism', fontsize=11)
ax.set_ylabel('Genus', fontsize=11)
ax.set_title('Genus x Resistance Mechanism Heatmap (Top 12 Genera)',
             fontsize=13, fontweight='bold', pad=12)
plt.xticks(rotation=35, ha='right', fontsize=9)
plt.yticks(fontsize=9)
fig.patch.set_facecolor('white')

plt.tight_layout()
plt.savefig('10_genus_mechanism_heatmap.png', dpi=180, bbox_inches='tight')
plt.show()
print("Plot 10 saved: 10_genus_mechanism_heatmap.png")


# ════════════════════════════════════════════════════════════
#  PLOT 10 — AMR Gene Diversity per Sample (Shannon Index)
# ════════════════════════════════════════════════════════════

def shannon_diversity(series):
    counts = series.value_counts()
    probs = counts / counts.sum()
    return -np.sum(probs * np.log(probs + 1e-10))

diversity = (filt.groupby('Sample')['AMR_Gene_Family']
             .apply(shannon_diversity).reindex(sample_order))

colors_div = plt.cm.RdYlGn(np.linspace(0.15, 0.85, len(diversity)))

fig, ax = plt.subplots(figsize=(12, 5))
bars = ax.bar(sample_order, diversity.values,
              color=colors_div, edgecolor='white', linewidth=0.8, width=0.65)

for bar, v in zip(bars, diversity.values):
    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
            f'{v:.2f}', ha='center', va='bottom',
            fontsize=8.5, fontweight='bold')

ax.set_xlabel('Sample', fontsize=11)
ax.set_ylabel('Shannon Diversity Index', fontsize=11)
ax.set_title('AMR Gene Family Diversity per Sample (Shannon Index)',
             fontsize=13, fontweight='bold', pad=12)
ax.set_facecolor('#f8f8f8')
fig.patch.set_facecolor('white')
ax.spines[['top', 'right']].set_visible(False)
ax.set_ylim(0, diversity.max() * 1.15)

plt.tight_layout()
plt.savefig('11_amr_diversity_shannon.png', dpi=180, bbox_inches='tight')
plt.show()
print("Plot 11 saved: 11_amr_diversity_shannon.png")


print("\nAll plots 8-11 generated successfully!")

