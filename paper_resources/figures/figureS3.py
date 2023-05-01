import json
import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator

THIS_SCRIPT_PATH = os.path.abspath(__file__)
DATA_PATH = os.path.join(os.path.dirname(THIS_SCRIPT_PATH), "data")
OUTPUT_FOLDER = os.path.join(os.path.dirname(THIS_SCRIPT_PATH), "output", "figureS3")


def main():
    PAE_TH = 50
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)

    pairs_validation_v2 = json.load(open(os.path.join(DATA_PATH, "benchmark1", "pairs_validation.json"), "r"))
    pairs_validation_v3 = json.load(open(os.path.join(DATA_PATH, "benchmark2", "pairs_validation.json"), "r"))

    v2_scores = []
    for jobname in pairs_validation_v2:
        for pair_scores in pairs_validation_v2[jobname]["pairs_scores"].values():
            scores_to_use = [i for i in pair_scores if i["pae_score"] > PAE_TH]
            if len(scores_to_use) == 0:
                continue
            v2_scores.append(max([i["dockq"] for i in scores_to_use]))

    v3_scores = []
    for jobname in pairs_validation_v3:
        for pair_scores in pairs_validation_v3[jobname]["pairs_scores"].values():
            scores_to_use = [i for i in pair_scores if i["pae_score"] > PAE_TH]
            if len(scores_to_use) == 0:
                continue
            v3_scores.append(max([i["dockq"] for i in scores_to_use]))

    print("medians", np.median(v2_scores), np.median(v3_scores))
    for DOCKQ_TH in [0.23, 0.49, 0.8]:
        print("dockq th:", DOCKQ_TH)
        print("v2", len(v2_scores), len([i for i in v2_scores if i > DOCKQ_TH]))
        print("v3", len(v3_scores), len([i for i in v3_scores if i > DOCKQ_TH]))

    data_to_plot = [v2_scores, v3_scores]

    fig, ax = plt.subplots(figsize=(2.25, 1.75), dpi=600)
    ax.xaxis.set_major_locator(MaxNLocator(6))
    ax.yaxis.set_major_locator(MaxNLocator(6))
    plt.yticks(fontsize=8)
    plt.xticks(fontsize=8)

    ax.violinplot(data_to_plot, showmedians=True)

    ax.set_xticks([1, 2])
    ax.set_xticklabels(["AFMv2\n(Benchmark 1)", "AFMv3\n(Benchmark 2)"])
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.ylabel('DockQ', fontsize=11)

    plt.savefig(os.path.join(OUTPUT_FOLDER, "FigS3.png"), bbox_inches='tight', dpi=300)


if __name__ == "__main__":
    main()
