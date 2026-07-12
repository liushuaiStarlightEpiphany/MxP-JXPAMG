import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Read data
data = []
with open('/tmp/opencode/exp_data.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        data.append(row)

# Groups to plot (same as paper)
groups = ['diag', 'ps1', 'ps2', 'ps3', 'general']
eps = '1e-4'

# Group colors
colors = {'diag': '#1f77b4', 'ps1': '#ff7f0e', 'ps2': '#2ca02c', 'ps3': '#d62728', 'general': '#9467bd'}
labels = {'diag': 'diag',
          'ps1': 'ps1',
          'ps2': 'ps2',
          'ps3': 'ps3',
          'general': 'general'}

fig, ax = plt.subplots(figsize=(8, 4.5))

for g in groups:
    errs = sorted([float(r['relerr'])
                   for r in data
                   if r['group'] == g and r['eps'] == eps])
    xs = np.arange(1, len(errs) + 1) / len(errs) * 100
    ax.semilogy(xs, errs, color=colors[g], label=labels[g], linewidth=1.2)

ax.set_xlabel('Matrix index sorted by error (%)', fontsize=11)
ax.set_ylabel('Relative error', fontsize=11)
ax.legend(fontsize=10, ncol=5, loc='upper center',
          bbox_to_anchor=(0.5, 1.12))
ax.set_xlim(0, 100)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/tmp/opencode/error_lines.png', dpi=200)
print("Saved to /tmp/opencode/error_lines.png")
