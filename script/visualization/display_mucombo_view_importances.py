import h5py
import numpy as np
from matplotlib import pyplot as plt
from sklearn.metrics import f1_score
from sklearn.model_selection import train_test_split
from sklearn.tree import DecisionTreeClassifier
from sklearn.utils import resample
from multimodal.boosting.combo import MuComboClassifier


def load_and_prepare(hdf5_path, labels_to_keep=(1, 2), view_keyword="View", label_key="Labels"):
    """
    Loads views and labels, filters classes, and concatenates features.
    Returns X, y, view indices, and view names.
    """
    with h5py.File(hdf5_path, "r") as f:
        view_ids = [k for k in f.keys() if view_keyword in k]
        X_views = [f[k][:] for k in view_ids]
        view_names = [f[k].attrs['name'] for k in view_ids]
        y = f[label_key][:]

    mask = np.isin(y, labels_to_keep)
    X = np.concatenate([v[mask] for v in X_views], axis=1)
    y = y[mask]

    # Compute MuCombo view_indices
    view_indices = [0] + np.cumsum([v.shape[1] for v in X_views]).tolist()
    return X, y, view_indices, view_names


def plot_view_importances(importances, view_names=None, top_k=3):
    """
    Displays view importance
    """
    importances = np.array(importances).mean(axis=0)
    sorted_idx = np.argsort(importances)[::-1]
    importances = importances[sorted_idx]

    if view_names is None:
        view_names = [f"View {i}" for i in range(len(importances))]
    else:
        view_names = [f"{name}" for name in view_names]
        view_names = np.array(view_names)[sorted_idx]

    colors = ["black" if i < top_k else "lightgray" for i in range(len(importances))]

    # Change values to display only a top subset of views
    n_views_to_keep = len(importances)

    plt.figure(figsize=(10, 6))
    bars = plt.barh(range(len(importances[:n_views_to_keep])), importances[:n_views_to_keep], color=colors)
    plt.yticks(range(len(importances[:n_views_to_keep])), view_names[:n_views_to_keep])
    plt.gca().invert_yaxis()  # The most important ones at the top

    # Adding score on the Figure
    for i, bar in enumerate(bars):
        width = bar.get_width()
        plt.text(width + 0.005, bar.get_y() + bar.get_height()/2,
                 f"{width:.3f}", va='center', fontsize=12)

    plt.xlabel("Importance Score", fontsize=16)
    plt.title("MuCombo Top 7 Views' Importance (Softmax Normalization) on EMF_MULTI_R", fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=20)

    plt.tight_layout()
    plt.show()


def evaluate_and_plot(X, y, view_indices, view_names,
                      max_depth=10, n_estimators=140,
                      n_splits=10, test_size=0.25, random_state=42, tau=0.25):
    """
    Compute MuCombo test F1 score and plot view importances
    """
    scores, importances = [], []
    for i in range(n_splits):
        rs = random_state + i
        Xt, Xv, yt, yv = train_test_split(X, y, test_size=test_size,
                                          stratify=y, random_state=rs)
        clf = MuComboClassifier(
            base_estimator=DecisionTreeClassifier(max_depth=max_depth),
            random_state=rs, n_estimators=n_estimators
        )
        clf.fit(Xt, yt, views_ind=view_indices)

        # F1 score
        y_pred = clf.predict(Xv)
        scores.append(f1_score(yv, y_pred, average="binary", pos_label=2))

        # Softmax with temperature
        w = clf.estimator_weights_alpha_.mean(axis=0)
        w = np.exp(w / tau)
        importances.append(w / w.sum())

    # Bootstrap confidence intervals 95%
    boot_means = [np.mean(resample(scores, replace=True, random_state=random_state))
                  for _ in range(1000)]
    lower, upper = np.percentile(boot_means, [2.5, 97.5])
    mean_f1, delta = np.mean(boot_means), (upper - lower) / 2
    print(f"F1 score (EMF): {mean_f1:.4f} ± {delta:.4f}")

    # Plotting importances
    plot_view_importances(importances, view_names)


def main():
    # Change the path with the right one if needed
    path = "../../data/dataset_compilation_to_hdf5/MuPPI_2025_14view_EMF.hdf5"
    X, y, view_indices, view_names = load_and_prepare(path)

    # Some stats
    print(f"X shape: {X.shape}")
    for cls, cnt in zip(*np.unique(y, return_counts=True)):
        print(f"Classe {cls}: {cnt / len(y):.2%}")
    evaluate_and_plot(X, y, view_indices, view_names)


if __name__ == '__main__':
    main()
