"""
This script makes a zoom on the five best classifiers with the highest score
"""

import pandas as pd
import matplotlib.pyplot as plt

# Load the file (replace the file name by the right one that you can find in your SUMMIT result directory)
df = pd.read_csv(
    "../path/to/your/summit/results/MuPPI_2025_14view_EMF_temp_filter-mean_on_10_iter-f1_score_p.csv",
    index_col=0)

# Keep only the 5 last ones (the best)
df = df.iloc[:, -5:]

# Extraction
train_scores = df.loc["Train"]
train_std = df.loc["Train STD"]
test_scores = df.loc["Test"]
test_std = df.loc["Test STD"]

# filtering the classifiers names
classifiers = df.columns
classifiers = [name.replace("-PPInetwork_topology", "-PPITopo") for name in classifiers]
classifiers = [name.replace("random_forest", "RF") for name in classifiers]
classifiers = [name.replace("weighted_linear_", "") for name in classifiers]
classifiers = [name.replace("imbalance_bagging", "IB") for name in classifiers]

x = range(len(classifiers))
bar_width = 0.35

plt.figure(figsize=(10, 6))

# Bars
bars_train = plt.bar([i - bar_width/2 for i in x], train_scores, width=bar_width, yerr=train_std,
                     alpha=0.7, label='Train', color="gray", capsize=5)
bars_test = plt.bar([i + bar_width/2 for i in x], test_scores, width=bar_width, yerr=test_std,
                    alpha=0.7, label='Test', color="black", capsize=5)

# Text under each bar
for i in x:
    plt.text(i - bar_width/2, -0.01, f"{train_scores[i]:.2f}±{train_std[i]:.2f}",
             ha='center', va='top', fontsize=10, rotation=0, color='gray')
    plt.text(i + bar_width/2, -0.01, f"{test_scores[i]:.2f}±{test_std[i]:.2f}",
             ha='center', va='top', fontsize=10, rotation=0, color='black')


plt.xticks(x, classifiers, rotation=45, fontsize=20, color='black')
plt.ylabel("F1 Score", fontsize=16, color='gray')
plt.title("Top Five Classifiers' Average F1 Scores on the EMF_MULTI_R task - SuMMIT Benchmark",
          fontsize=16, color='black')
plt.legend(fontsize=14)
plt.ylim(bottom=-0.1)  # to leave space for the text under the bars
plt.tight_layout()
plt.show()
