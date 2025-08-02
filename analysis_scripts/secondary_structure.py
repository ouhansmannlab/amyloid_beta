#This is the final code that used to plot the figures in the manuscript
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

# Load trajectory
traj = md.load('recheck/rep1_ss.xtc', top='rep1_ss.gro')

# Compute DSSP with simplified secondary structures
dssp = md.compute_dssp(traj, simplified=True)

# Define DSSP codes
dssp_codes = {
    'H': 0,  # Helix
    'E': 1,  # Strand
    'C': 2,  # Coil
}

# Convert DSSP to numeric array
dssp_numeric = np.array([[dssp_codes[code] for code in frame] for frame in dssp])

# Calculate residue-wise percentage of secondary structures
num_frames = dssp_numeric.shape[0]
num_residues = dssp_numeric.shape[1]
percentages = np.zeros((len(dssp_codes), num_residues))

for i in range(num_residues):
    residue_dssp = dssp_numeric[:, i]
    for code, numeric in dssp_codes.items():
        percentages[numeric, i] = np.sum(residue_dssp == numeric) / num_frames * 100

# Prepare data for CSV export
structure_labels = np.array(list(dssp_codes.keys())).reshape(-1, 1)
data = np.hstack((structure_labels, percentages))

# Export to CSV
header = 'Structure,' + ','.join([f'Residue_{i+1}' for i in range(num_residues)])
np.savetxt('residue_wise_percentage_secondary_structure.csv', data, delimiter=',', fmt='%s', header=header, comments='')

# Plotting secondary structure over time
colors = ['red', 'blue', 'springgreen']
cmap = ListedColormap(colors[:len(dssp_codes)])
boundaries = np.arange(len(dssp_codes) + 1)  # [0, 1, 2, 3]
norm = BoundaryNorm(boundaries, cmap.N)

plt.figure(dpi=300)
img = plt.imshow(dssp_numeric.T, aspect='auto', cmap=cmap, norm=norm, interpolation='nearest')

# Add colorbar with centered ticks
cbar = plt.colorbar(img, ticks=np.arange(len(dssp_codes)) + 0.5)
cbar.set_ticklabels(list(dssp_codes.keys()))
cbar.ax.tick_params(labelsize=18)
cbar.set_label('Secondary Structure', fontsize=18)

# Axis labels and ticks
plt.xlabel('Time (ns)', fontsize=18)
plt.ylabel('Residue Index', fontsize=18)
plt.title('RC1', fontsize=20)
plt.xticks([1000, 2000, 3000])

residue_indices = np.arange(num_residues)
residue_labels = [str(i+1) if ((i+1) % 10 == 0 or (i+1) == 1) else '' for i in range(num_residues)]
plt.yticks(residue_indices, residue_labels, fontsize=18)
plt.xticks(fontsize=18)

plt.grid(False)
plt.gca().set_facecolor('white')
plt.gca().invert_yaxis()
plt.savefig('manuscript_fig/rep1_secondary_structure.png', bbox_inches='tight', facecolor='white')
plt.show()

# Plotting residue-wise percentage of secondary structures
fig, ax = plt.subplots(dpi=300)
width = 0.8
x = np.arange(num_residues)

for i, (structure, color) in enumerate(zip(dssp_codes.keys(), colors)):
    ax.bar(x, percentages[i], width, label=structure, color=color, bottom=np.sum(percentages[:i], axis=0))

ax.set_xlabel('Residue Index', fontsize=16)
ax.set_ylabel('Percentage', fontsize=16)
ax.set_title('Residue-wise Percentage of Secondary Structure', fontsize=16)
ax.set_xticks([0, 9, 19, 29, 39])
ax.set_xticklabels(['1', '10', '20', '30', '40'], fontsize=16)
ax.tick_params(axis='y', labelsize=16)
ax.legend(fontsize=14)

plt.grid(False)
plt.gca().set_facecolor('white')
plt.savefig('manuscript_fig/rep1_residue_wise_percentage_secondary_structure.png', bbox_inches='tight', facecolor='white')
plt.show()

