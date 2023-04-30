import json
import os
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator


THIS_SCRIPT_PATH = os.path.abspath(__file__)
DATA_PATH = os.path.join(os.path.dirname(THIS_SCRIPT_PATH), "data")
OUTPUT_FOLDER = os.path.join(os.path.dirname(THIS_SCRIPT_PATH), "output", "figure2")


def draw_scatter(name1, name2, x, y, output_path):
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

    equation = f'$\\rho = {np.corrcoef(x, y)[0][1]:.2f}$'
    plt.annotate(equation, xy=(0.05, 0.9), xycoords='axes fraction', fontsize=8,
                 bbox=dict(boxstyle='square', facecolor='white', edgecolor="black", lw=0.5, pad=0.2))

    plt.savefig(output_path, bbox_inches='tight', dpi=300)


def get_connectivity_groups(pairs_scores, all_names, dockq_th: float = 0.23, pae_th: float = 0):
    subunits_groups = {subunit_name: i for i, subunit_name in enumerate(all_names)}
    for subunit_pair, scores in pairs_scores.items():
        subunit_name1, subunit_name2 = subunit_pair.split("_plus_")
        if subunits_groups[subunit_name1] == subunits_groups[subunit_name2]:
            continue
        for score in scores:
            if score["dockq"] > dockq_th and score["pae_score"] > pae_th:
                group_to_change = subunits_groups[subunit_name2]
                for subunit_name in subunits_groups.keys():
                    if subunits_groups[subunit_name] == group_to_change:
                        subunits_groups[subunit_name] = subunits_groups[subunit_name1]
                break
    group_to_subunit = defaultdict(list)
    for subunit, group_id in subunits_groups.items():
        group_to_subunit[group_id].append(subunit)
    connectivity_groups = sorted(list(group_to_subunit.values()), key=len, reverse=True)
    return connectivity_groups


def main():
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)

    # load benchmark info
    pdb_to_subunits = {}
    benchmark_path = os.path.join(DATA_PATH, "benchmark1")

    for filename in os.listdir(os.path.join(benchmark_path, "pdb_infos")):
        jobname = filename.split(".")[0]
        if not filename.endswith(".json"):
            continue
        pdb_to_subunits[jobname] = json.load(open(os.path.join(benchmark_path, "pdb_infos", filename), "rb"))

    print("Benchmark size ", len(pdb_to_subunits), "PDB ids: ", list(pdb_to_subunits.keys()))

    # load combfold results
    combfold_results = json.load(open(os.path.join(benchmark_path, "combfold_results.json"), "r"))
    combfold_results = {k: v for k, v in combfold_results.items() if v is not None and k in pdb_to_subunits}
    combfold_parsed_results = {}
    for jobname, results in combfold_results.items():
        if not results["scores"]:
            continue
        tm_and_scores = [(result['weighted_trans_score'], result['tm_score']) for result in results["scores"].values()]
        combfold_parsed_results[jobname] = [i[1] for i in sorted(tm_and_scores, reverse=True)]

    pairs_validation = json.load(open(os.path.join(benchmark_path, "pairs_validation.json"), "r"))

    print("CombFold able to assemble ", len(combfold_parsed_results))
    print("CombFold able to assemble correctly (TM-score > 0.7) ", len([i for i in combfold_parsed_results.values()
                                                                        if max(i) > 0.7]))

    # load AFMv2 results
    af2_results_path = os.path.join(benchmark_path, "AFMv2_results.json")
    af2_results = {}

    validation_results = json.load(open(af2_results_path, "r"))
    for jobname in validation_results:
        if jobname not in pdb_to_subunits:
            continue
        af2_results[jobname] = []
        for pdb_name in sorted(validation_results[jobname].keys()):
            af2_results[jobname].append(validation_results[jobname][pdb_name]["tm_score"])
    print("AFMv2 able to assemble ", len([i for i in af2_results.values() if i]))

    # Fig 2A - bar graph - Combfold vs. AFMv2
    name1, name2 = "CombFold", "AFMv2"
    th_high, th_accept = 0.8, 0.7
    bar_width = 0.3
    fig, ax = plt.subplots()

    labels = []
    for count, max_t in enumerate([1, 5, 10]):
        labels.append(f"Top {max_t}")
        results1 = [max(v[:max_t], default=0) for v in combfold_parsed_results.values()]
        results2 = [max(v[:max_t], default=0) for v in af2_results.values()]
        res1_high = len([i for i in results1 if i >= th_high])
        res1_acceptable = len([i for i in results1 if th_high > i >= th_accept])
        res2_high = len([i for i in results2 if i >= th_high])
        res2_acceptable = len([i for i in results2 if th_high > i >= th_accept])

        print(name1, max_t, res1_high, res1_acceptable)
        print(name2, max_t, res2_high, res2_acceptable)

        # ax.bar(count, [res1_high / len(pdb_to_subunits)], color='#1f78b4', width=bar_width, edgecolor='grey',
        #        label=f"{name1}\nHigh")
        # ax.bar(count, [res1_acceptable / len(pdb_to_subunits)], color='#1f78b4', alpha=0.5, width=bar_width,
        #        edgecolor='grey', bottom=[res1_high / len(pdb_to_subunits)], label=f"{name1}\nAcceptable")
        # ax.bar(count + bar_width, [res2_high / len(pdb_to_subunits)], color='#ff7f00', width=bar_width, edgecolor='grey',
        #        label=f"{name2}\nHigh")
        # ax.bar(count + bar_width, [res2_acceptable / len(pdb_to_subunits)], color='#ff7f00', alpha=0.5, width=bar_width,
        #        edgecolor='grey', bottom=[res2_high / len(pdb_to_subunits)], label=f"{name2}\nAcceptable")
        ax.bar(count, [res1_high / len(pdb_to_subunits)], color='#1f78b4', width=bar_width, edgecolor='grey',
               label=f"{name1} High")
        ax.bar(count, [res1_acceptable / len(pdb_to_subunits)], color='#1f78b4', alpha=0.5, width=bar_width,
               edgecolor='grey', bottom=[res1_high / len(pdb_to_subunits)], label=f"{name1} Acceptable")
        ax.bar(count + bar_width, [res2_high / len(pdb_to_subunits)], color='#ff7f00', width=bar_width,
               edgecolor='grey',
               label=f"{name2} High")
        ax.bar(count + bar_width, [res2_acceptable / len(pdb_to_subunits)], color='#ff7f00', alpha=0.5, width=bar_width,
               edgecolor='grey', bottom=[res2_high / len(pdb_to_subunits)], label=f"{name2} Acceptable")

    plt.ylabel('Success rate', fontsize=11)
    plt.xticks([i + bar_width / 2 for i in range(len(labels))], labels, fontsize=8)
    plt.yticks(np.arange(0.1, 0.9, 0.1), fontsize=8)

    handles, labels = plt.gca().get_legend_handles_labels()
    new_labels, new_handles = [], []
    for handle, label in zip(handles, labels):
        if label not in new_labels:
            new_labels.append(label)
            new_handles.append(handle)
    ax.legend(new_handles, new_labels, bbox_to_anchor=(-0.35, -0.7), loc='lower left', fontsize=8, ncol=2,
              handletextpad=0.2, columnspacing=0.25)
    # ax.legend(new_handles, new_labels, bbox_to_anchor=(-0.4, 1.05), loc='lower left',  fontsize=6, ncol=2)
    fig.set_size_inches(2, 1.5)
    fig.set_dpi(300)

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)

    plt.savefig(os.path.join(OUTPUT_FOLDER, "Fig2A.png"), bbox_inches='tight', dpi=300)

    # Fig 2B - Confidence - TM-score
    confidence, tm_results = [], []

    for jobname, results in combfold_results.items():
        if len(results["scores"]) == 0:
            continue
        conf_key = "weighted_trans_score"
        best_confidence = max(list(results["scores"].values()), key=lambda x: x[conf_key])
        confidence.append(best_confidence[conf_key])
        tm_results.append(best_confidence["tm_score"])
    print("confidence correlation:", np.corrcoef(confidence, tm_results)[0][1])
    draw_scatter("Predicted confidence", "TM-score", confidence, tm_results, os.path.join(OUTPUT_FOLDER, "Fig2B.png"))

    # Fig 2C - Weighted connectivity Ratio - TM-score
    tm_results = []
    pairwise_connectivity = []
    MAX_T = 10

    for jobname in sorted(pdb_to_subunits.keys()):
        if jobname not in pairs_validation:
            continue
        tm_results.append(max(combfold_parsed_results.get(jobname, [0])[:MAX_T]))

        connectivity_groups = get_connectivity_groups(pairs_validation[jobname]["pairs_scores"],
                                                      list(pairs_validation[jobname]["domains_scores"].keys()),
                                                      dockq_th=0.23, pae_th=0)
        total_count = connected_count = 0
        for domain_name in pairs_validation[jobname]["domains_scores"].keys():
            res_num = len(pdb_to_subunits[jobname][domain_name.split("_")[0]]["sequence"])
            total_count += res_num
            if domain_name in connectivity_groups[0] \
                    and pairs_validation[jobname]["domains_scores"][domain_name]["rmsd"] < 10:
                connected_count += res_num
        pairwise_connectivity.append(connected_count / total_count)
        if pairwise_connectivity[-1] >= 0.6 and tm_results[-1] < 0.7:
            print("****", jobname, pairwise_connectivity[-1], tm_results[-1])
    print(np.corrcoef(pairwise_connectivity, tm_results)[0][1])
    print("probably possible", len([i for i in pairwise_connectivity if i >= 0.6]))
    draw_scatter("Pairwise Connectivity", "TM-score", pairwise_connectivity, tm_results,
                 os.path.join(OUTPUT_FOLDER, "Fig2C.png"))
    print("Pairwise Connectivity correlation:", np.corrcoef(pairwise_connectivity, tm_results)[0][1])

    # Fig 2D - scatter plot - Combfold vs. AFMv2
    MAX_T = 1
    combfold_tm_results, afm2_tm_results = [], []
    for jobname in pdb_to_subunits:
        combfold_tm_results.append(max(combfold_parsed_results.get(jobname, [])[:MAX_T], default=0))
        afm2_tm_results.append(max(af2_results.get(jobname, [])[:MAX_T], default=0))
        if afm2_tm_results[-1] > 0.7 > combfold_tm_results[-1]:
            print(jobname, "better on AFMv2", afm2_tm_results[-1], combfold_tm_results[-1])

    fig, ax = plt.subplots(figsize=(2.25, 1.75), dpi=600)
    ax.scatter(combfold_tm_results, afm2_tm_results, alpha=0.4, s=10, edgecolor='none')
    ax.set_xlabel("CombFold TM-score", fontsize=11)
    ax.set_ylabel("AFMv2 TM-score", fontsize=11)
    plt.yticks(np.arange(0, 1.2, 0.2), fontsize=8)
    plt.xticks(np.arange(0, 1.2, 0.2), fontsize=8)

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.plot([0, 1], [0, 1], 'r--')

    plt.savefig(os.path.join(OUTPUT_FOLDER, "Fig2D.png"), bbox_inches='tight', dpi=300)


if __name__ == "__main__":
    main()