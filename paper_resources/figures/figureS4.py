import json
import math
import os

import matplotlib.pyplot as plt
import numpy as np
import scipy


THIS_SCRIPT_PATH = os.path.abspath(__file__)
DATA_PATH = os.path.join(os.path.dirname(THIS_SCRIPT_PATH), "data")
OUTPUT_FOLDER = os.path.join(os.path.dirname(THIS_SCRIPT_PATH), "output", "figureS4")


def main():
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)

    # load benchmark info
    pdb_to_subunits = {}
    benchmark_path = os.path.join(DATA_PATH, "benchmark2")

    for filename in os.listdir(os.path.join(benchmark_path, "pdb_infos")):
        jobname = filename.split(".")[0]
        if not filename.endswith(".json"):
            continue
        pdb_to_subunits[jobname] = json.load(open(os.path.join(benchmark_path, "pdb_infos", filename), "rb"))

    print("Benchmark size ", len(pdb_to_subunits), "PDB ids: ", list(pdb_to_subunits.keys()))

    # fig a - ICS scores
    afm_ics = {'7QRU': 0.8875607022786702, '7OZN': 0, '8BEL': 0.7880794701986755, '7ZKQ': 0.3887775551102204, '7E8T': 0, '7T3B': 0.6454709058188361, '7OBA': 0, '8ADN': 0, '7ARC': 0.4630650496141124, '7XHO': 0, '8CTE': 0, '8A3T': 0, '7T2R': 0.8136942675159236, '8E9G': 0, '8ADL': 0, '7UIC': 0, '8HIL': 0.663960306720794, '7PKN': 0.735445836403832, '7USE': 0.823474615742897, '8FEE': 0.7976314872866598, '7P2Y': 0, '7QVE': 0, '7WFF': 0, '8A5O': 0.7162534435261708, '8F5O': 0}
    combfold_ics = {'7QRU': 0.8370890284447728, '7OZN': 0, '8BEL': 0.6606606606606606, '7ZKQ': 0.1755062680810029, '7E8T': 0.26416168516937094, '7T3B': 0.659914712153518, '7OBA': 0.49898511502029774, '8ADN': 0.39416432045047345, '7ARC': 0.360079387581514, '7XHO': 0.12848370335380258, '8CTE': 0.31272210376687987, '8A3T': 0.19494769073743298, '7T2R': 0.7802621434078643, '8E9G': 0.5614627285513363, '8ADL': 0, '7UIC': 0.16527390900649955, '8HIL': 0.35519387266634755, '7PKN': 0.22074074074074074, '7USE': 0.5253557095950383, '8FEE': 0.06477409036385445, '7P2Y': 0.34454358270418667, '7QVE': 0.2577530719719134, '7WFF': 0.5515331355093966, '8A5O': 0.24546424759871932, '8F5O': 0.0031237797735259665}
    # relaxed:
    # combfold_ics = {'7QRU': 0.7505800464037122, '7OZN': 0, '8BEL': 0.6078199052132701, '7ZKQ': 0.17944250871080136, '7E8T': 0.2566658360328931, '7T3B': 0.6181506849315069, '7OBA': 0.47959327778562344, '8ADN': 0, '7ARC': 0.3612124261099233, '7XHO': 0, '8CTE': 0, '8A3T': 0, '7T2R': 0.7269395235186317, '8E9G': 0.5418730832743571, '8ADL': 0, '7UIC': 0.16054654141759178, '8HIL': 0.34356084324775626, '7PKN': 0.21782639050500346, '7USE': 0.49237472766884527, '8FEE': 0.08069277701239914, '7P2Y': 0.37913858775271025, '7QVE': 0, '7WFF': 0.557182820470627, '8A5O': 0.2362052274927396, '8F5O': 0}
    afm_results, combfold_results = [], []
    for pdb_id in afm_ics:
        afm_results.append(afm_ics[pdb_id])
        combfold_results.append(combfold_ics[pdb_id])

    fig, ax = plt.subplots(figsize=(2.25, 1.75), dpi=600)
    ax.scatter(combfold_results, afm_results, alpha=0.4, s=10, edgecolor='none')
    ax.set_xlabel("CombFold ICS", fontsize=11)
    ax.set_ylabel("AFMv3 ICS", fontsize=11)
    plt.yticks(np.arange(0, 1.2, 0.2), fontsize=8)
    plt.xticks(np.arange(0, 1.2, 0.2), fontsize=8)

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.plot([0, 1], [0, 1], 'r--')

    plt.savefig(os.path.join(OUTPUT_FOLDER, "FigS4a.png"), bbox_inches='tight', dpi=300)

    # Fig b - Prodigy scores
    prodigy_scores = json.load(open(os.path.join(benchmark_path, "prodigy_results.json"), "r"))

    joined_scores_targets, joined_scores_models = [], []
    for pdb_id, scores in prodigy_scores.items():
        if scores is None or len(scores) == 0:
            continue
        for score_key, score in scores.items():
            if score[1] is None:
                continue

            if min(score[0], score[1]) < 1e-30:
                print("Skipping ", score[0], score[1], "for ", pdb_id, "as it is too small")
                continue
            if abs(math.log(score[0], 10) - math.log(score[1], 10)) > 10:
                print("Large diff ", pdb_id, score_key, score[0], score[1])
            joined_scores_targets.append(score[0])
            joined_scores_models.append(score[1])

        print("values", len(scores), joined_scores_targets[-1 * len(scores):], joined_scores_models[-1 * len(scores):])
        print("correlation for ", pdb_id, np.corrcoef(joined_scores_targets[-1 * len(scores):],
                                                      joined_scores_models[-1 * len(scores):])[0][1])

    fig, ax = plt.subplots(figsize=(2.25, 1.75), dpi=600)
    ax.scatter(joined_scores_targets, joined_scores_models, alpha=0.4, s=10, edgecolor='none')
    ax.set_xlabel("Experimental structure", fontsize=11)
    ax.set_ylabel("CombFold Top-1 model", fontsize=11)

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.plot([0, 1], [0, 1], 'r--', lw=0.5)

    plt.xscale('log')
    plt.yscale('log')

    plt.xlim(10 ** -30, 1)
    plt.ylim(10 ** -30, 1)

    corr = scipy.stats.spearmanr(joined_scores_targets, joined_scores_models)[0]

    equation = f'$\\rho = {corr:.2f}$'

    plt.annotate(equation, xy=(0.05, 0.9), xycoords='axes fraction', fontsize=8,
                 bbox=dict(boxstyle='square', facecolor='white', edgecolor="black", lw=0.5, pad=0.2))

    plt.savefig(os.path.join(OUTPUT_FOLDER, "FigS4b.png"), bbox_inches='tight', dpi=300)

    # clashes:
    not_relaxed, relaxed = [], []
    for line in open(os.path.join(benchmark_path, "clashscores.csv"), "r").read().split("\n")[1:]:
        if not line:
            continue
        pdb_id, score_not_relaxed, score_relaxed = line.split(",")[:3]
        not_relaxed.append(float(score_not_relaxed))
        relaxed.append(float(score_relaxed))

    data_to_plot = [not_relaxed, relaxed]

    fig, ax = plt.subplots(figsize=(2.25, 1.75), dpi=600)
    plt.yticks(fontsize=11)
    plt.xticks(fontsize=11)

    ax.violinplot(data_to_plot, showmedians=True)

    ax.set_xticks([1, 2])
    ax.set_xticklabels(["No relaxation", "Relaxed"])
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.ylabel('MolProbity clashscore', fontsize=11)

    plt.savefig(os.path.join(OUTPUT_FOLDER, "FigS4c.png"), bbox_inches='tight', dpi=300)


if __name__ == "__main__":
    main()