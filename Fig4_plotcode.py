import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ranksums

def create_boxplot(data, title, xlabel, ylabel, xticklabels, colormap, colormap_range, figure_size=(8, 6), font_size=16, line_width=2):
    fig, ax = plt.subplots(figsize=figure_size)
    ax.set_xlabel(xlabel, fontsize=font_size)
    ax.set_ylabel(ylabel, fontsize=font_size)
    ax.set_title(title, fontsize=font_size + 2)

    cmap = plt.get_cmap(colormap)
    colors = cmap(np.linspace(colormap_range[0], colormap_range[1], len(data)))

    bplot = plt.boxplot(data, positions=range(len(data)), patch_artist=True)

    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_edgecolor('black')
        patch.set_linewidth(line_width)
    for whisker in bplot['whiskers']:
        whisker.set(color='black', linewidth=line_width)
    for cap in bplot['caps']:
        cap.set(color='black', linewidth=line_width)
    for median in bplot['medians']:
        median.set(color='yellow', linewidth=line_width)
    for flier in bplot['fliers']:
        flier.set(marker='o', color='red', alpha=0.5, markersize=6)

    ax.grid(True, linestyle='--', linewidth=0.6, alpha=0.7)
    ax.set_xticks(range(len(data)))
    ax.set_xticklabels(xticklabels, fontsize=font_size)
    ax.tick_params(axis='both', which='major', labelsize=font_size)

    plt.tight_layout()
    plt.show()

# Volume of Area Cleaned Through Activation (Using Gray Color Scheme)
cleaning_bis = [[43.87,41.86,21.44,44.55,7.35,26.81,46.30],[50.71,45.19,54.85,60.74,20.65,36.73,60.74],[59.15,46.81,79,73.80,27.97,94.93,75],[66.89,54.58,100,90.34,34.85,100,81]]
create_boxplot(
    cleaning_bis,
    title=" ",
    xlabel=" ",
    ylabel="Cleaned Area (%)",
    xticklabels=['Burst 1', 'Burst 2', 'Burst 3', 'Burst 4'],
    colormap='Greys',
    colormap_range=(1, 0.3),
    font_size=26,
    line_width=2
)


# Data (replace with different data each time)
sample1 = [50.71,45.19,54.85,60.74,20.65,36.73,60.74]
sample2 = [59.15,46.81,79,73.80,27.97,94.93,75]

# Perform the test,
#‘less’: the distribution underlying sample1 is stochastically less than the distribution underlying sample2
#‘greater’: the distribution underlying sample1 is stochastically greater than the distribution underlying sample2.
stat, p_value = ranksums(sample1, sample2, alternative ='less')

print("Test statistic:", stat)
print("p-value:", p_value)



