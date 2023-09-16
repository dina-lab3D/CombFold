import json
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator

THIS_SCRIPT_PATH = os.path.abspath(__file__)
DATA_PATH = os.path.join(os.path.dirname(THIS_SCRIPT_PATH), "data")
OUTPUT_FOLDER = os.path.join(os.path.dirname(THIS_SCRIPT_PATH), "output", "figureS6")


def score_to_pae(score):
    # score = 100 - (1 / 4) * (full_pae_avg ** 2)
    return np.sqrt((100 - score) * 4)


def to_pae_scores(scores):
    return [score_to_pae(score) for score in scores]


def draw_scatter(name1, name2, x, y, output_path, y_lim=None):
    fig, ax = plt.subplots(figsize=(2.25, 1.75), dpi=600)
    ax.scatter(x, y, alpha=0.4, s=10, edgecolor='none', color="#1f78b4")
    ax.xaxis.set_major_locator(MaxNLocator(6))
    ax.yaxis.set_major_locator(MaxNLocator(6))
    plt.yticks(fontsize=8)
    plt.xticks(fontsize=8)
    plt.xlabel(name1, fontsize=11)
    plt.ylabel(name2, fontsize=11)
    p = np.poly1d(np.polyfit(x, y, 1))
    x_lin_space = np.linspace(min(x), max(x), 100)
    plt.plot(x_lin_space, p(x_lin_space), color="#ff7f00", linewidth=1, linestyle="--")
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)

    if y_lim is not None:
        plt.ylim(y_lim)

    equation = f'$\\rho = {np.corrcoef(x, y)[0][1]:.2f}$'
    # plt.annotate(equation, xy=(0.05, 0.9), xycoords='axes fraction', fontsize=8,
    #              bbox=dict(boxstyle='square', facecolor='white', edgecolor="black", lw=0.5, pad=0.2))
    plt.annotate(equation, xy=(0.7, 0.9), xycoords='axes fraction', fontsize=8,
                 bbox=dict(boxstyle='square', facecolor='white', edgecolor="black", lw=0.5, pad=0.2))

    plt.savefig(output_path, bbox_inches='tight', dpi=300)


def main():
    # load benchmark info
    benchmark_path = os.path.join(DATA_PATH, "benchmark2")
    scores_by_pdb = json.load(open(os.path.join(benchmark_path, "rmsd_validate.json"), "r"))
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)

    rmsds = []
    iptms = []
    iplddt = []
    weighted_iplddt = []
    full_pae = []
    processed_pae = []
    """
    {'pdb': 'HDF_NDF_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_000.pdb', 'iptm': 0.55, 
    'rmsd': 9.05883774873516, 's1_plddt': 76.62620098039216, 's1_iplddt': 76.20222222222222, 
    's1_pae': 10.695484849577085, 's2_plddt': 89.00289827255278, 's2_iplddt': 74.16875, 's2_pae': 5.845915539656869, 
    'full_pae': 17.015451902450042, 'interface_pae': 8.971527777777778, 'interface1_size': 9, 'interface2_size': 8} 
    """
    for pdb_id, scores in scores_by_pdb.items():
        for score in scores:
            rmsds.append(score["rmsd"])
            iptms.append(score["iptm"])
            iplddt.append((score["s1_iplddt"] + score["s2_iplddt"]) / 2)
            weighted_iplddt.append((score["interface1_size"] * score["s1_iplddt"]
                                    + score["interface2_size"] * score["s2_iplddt"])
                                   / (score["interface1_size"] + score["interface2_size"]))
            processed_pae.append(max([1, 100 - (1 / 4) * (score["full_pae"] ** 2)]))
            full_pae.append(score["full_pae"])

    draw_scatter("RMSD", "ipTM", rmsds, iptms, os.path.join(OUTPUT_FOLDER, "FigS6a.png"), y_lim=(0, 1))
    draw_scatter("RMSD", "ipLDDT", rmsds, iplddt, os.path.join(OUTPUT_FOLDER, "FigS6b.png"), y_lim=(0, 100))
    draw_scatter("RMSD", "Normalized ipLDDT", rmsds, weighted_iplddt, os.path.join(OUTPUT_FOLDER, "FigS6c.png"),
                 y_lim=(0, 100))
    draw_scatter("RMSD", "Full PAE", rmsds, full_pae, os.path.join(OUTPUT_FOLDER, "FigS6d.png"), y_lim=(0, 35))
    draw_scatter("RMSD", "CombFold Score", rmsds, processed_pae, os.path.join(OUTPUT_FOLDER, "FigS6e.png"),
                 y_lim=(0, 100))

    # Fig F
    all_pae_scores = json.load(open(os.path.join(DATA_PATH, "benchmark2", "pae_data", "all_pae_scores.json"), "r"))
    used_pae_scores = json.load(open(os.path.join(DATA_PATH, "benchmark2", "pae_data", "used_pae_scores.json"), "r"))
    top1_used_pae_scores = json.load(open(os.path.join(DATA_PATH, "benchmark2", "pae_data",
                                                       "top1_used_pae_scores.json"), "r"))

    data_to_plot = [to_pae_scores(all_pae_scores), to_pae_scores(top1_used_pae_scores)]

    fig, ax = plt.subplots(figsize=(2.25, 1.75), dpi=600)
    ax.xaxis.set_major_locator(MaxNLocator(6))
    ax.yaxis.set_major_locator(MaxNLocator(6))
    plt.yticks(fontsize=8)
    plt.xticks(fontsize=8)

    ax.violinplot(data_to_plot, showmedians=True)

    ax.set_xticks([1, 2])
    ax.set_xticklabels(["All Generated", "Used in Top1 Results"])
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.ylabel('PAE scores', fontsize=11)

    print("medians: ", np.median(all_pae_scores), np.median(top1_used_pae_scores))
    print("means: ", np.mean(all_pae_scores), np.mean(top1_used_pae_scores))

    plt.savefig(os.path.join(OUTPUT_FOLDER, "FigS6f.png"), bbox_inches='tight', dpi=300)


if __name__ == "__main__":
    main()