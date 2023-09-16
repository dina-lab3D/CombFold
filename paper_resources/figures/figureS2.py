import json
import os

import matplotlib.pyplot as plt
import numpy as np


THIS_SCRIPT_PATH = os.path.abspath(__file__)
DATA_PATH = os.path.join(os.path.dirname(THIS_SCRIPT_PATH), "data")
OUTPUT_FOLDER = os.path.join(os.path.dirname(THIS_SCRIPT_PATH), "output", "figureS2")


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

    # load combfold results
    combfold_results = json.load(open(os.path.join(benchmark_path, "combfold_results.json"), "r"))
    combfold_results = {k: v for k, v in combfold_results.items() if v is not None and k in pdb_to_subunits}
    combfold_parsed_results = {}
    for jobname, results in combfold_results.items():
        if not results["scores"]:
            continue
        tm_and_scores = [(result['weighted_trans_score'], result['tm_score']) for result in results["scores"].values()]
        combfold_parsed_results[jobname] = [i[1] for i in sorted(tm_and_scores, reverse=True)]

    print("CombFold able to assemble ", len(combfold_parsed_results))
    print("CombFold able to assemble correctly (TM-score > 0.7) ", len([i for i in combfold_parsed_results.values()
                                                                        if max(i) > 0.7]))

    # load AFMv3 results
    afm3_results_path = os.path.join(benchmark_path, "AFMv3_results.json")
    afm3_results = {}

    validation_results = json.load(open(afm3_results_path, "r"))
    all_jobnames = {i.split("_")[0] for i in validation_results}

    for jobname in all_jobnames:
        all_jobs = [i for i in validation_results if i.startswith(jobname)]
        all_jobs_ranked = [i for i in all_jobs if "rank" in i]
        if all_jobs_ranked:
            all_jobs = all_jobs_ranked
        all_jobs.sort()
        afm3_results[jobname] = [validation_results[filename]["tm_score"] for filename in all_jobs]

    print("AFMv3 able to assemble ", len([i for i in afm3_results.values() if i]))

    # load rossetta results
    rosetta_results_path = os.path.join(benchmark_path, "rosetta_results.json")
    rosetta_results_json = json.load(open(rosetta_results_path, "r"))
    rosetta_results = {jobname: [i["tm_score"] for i in results.values()] if results else []
                       for jobname, results in rosetta_results_json.items()}

    print("Rosetta able to assemble ", len([i for i in rosetta_results.values() if i]))

    # Create bar plot
    name1, name2, name3 = "CombFold", "AFMv3", "RosettaFold2"
    th_high, th_accept = 0.8, 0.7
    bar_width = 0.3
    fig, ax = plt.subplots()

    labels = []
    for count, max_t in enumerate([1, 5, 10]):
        labels.append(f"Top {max_t}")
        results1 = [max(v[:max_t], default=0) for v in combfold_parsed_results.values()]
        results2 = [max(v[:max_t], default=0) for v in afm3_results.values()]
        results3 = [max(v[:max_t], default=0) for v in rosetta_results.values()]
        res1_high = len([i for i in results1 if i >= th_high])
        res1_acceptable = len([i for i in results1 if th_high > i >= th_accept])
        res2_high = len([i for i in results2 if i >= th_high])
        res2_acceptable = len([i for i in results2 if th_high > i >= th_accept])
        res3_high = len([i for i in results3 if i >= th_high])
        res3_acceptable = len([i for i in results3 if th_high > i >= th_accept])

        print(name1, max_t, res1_high, res1_acceptable)
        print(name2, max_t, res2_high, res2_acceptable)
        print(name3, max_t, res3_high, res3_acceptable)

        ax.bar(count, [res1_high / len(pdb_to_subunits)], color='#1f78b4', width=bar_width, edgecolor='grey',
               label=f"{name1}\nHigh")
        ax.bar(count, [res1_acceptable / len(pdb_to_subunits)], color='#1f78b4', alpha=0.5, width=bar_width,
               edgecolor='grey', bottom=[res1_high / len(pdb_to_subunits)], label=f"{name1}\nAcceptable")
        ax.bar(count + bar_width, [res2_high / len(pdb_to_subunits)], color='#ff7f00', width=bar_width, edgecolor='grey',
               label=f"{name2}\nHigh")
        ax.bar(count + bar_width, [res2_acceptable / len(pdb_to_subunits)], color='#ff7f00', alpha=0.5, width=bar_width,
               edgecolor='grey', bottom=[res2_high / len(pdb_to_subunits)], label=f"{name2}\nAcceptable")
        ax.bar(count + 2 * bar_width, [res3_high / len(pdb_to_subunits)], color='#33a02c', width=bar_width, edgecolor='grey',
               label=f"{name3}\nHigh")
        ax.bar(count + 2 * bar_width, [res3_acceptable / len(pdb_to_subunits)], color='#33a02c', alpha=0.5, width=bar_width,
               edgecolor='grey', bottom=[res3_high / len(pdb_to_subunits)], label=f"{name3}\nAcceptable")

    plt.ylabel('Success rate', fontsize=11)
    plt.xticks([i + bar_width / 2 for i in range(len(labels))], labels, fontsize=8)
    plt.yticks(np.arange(0.1, 0.9, 0.1), fontsize=8)

    handles, labels = plt.gca().get_legend_handles_labels()
    new_labels, new_handles = [], []
    for handle, label in zip(handles, labels):
        if "High" not in label:
            continue
        small_label = label.replace("\nHigh", "")
        if small_label not in new_labels:
            new_labels.append(small_label)
            new_handles.append(handle)
    ax.legend(new_handles, new_labels, bbox_to_anchor=(1, 0), loc='lower left',  fontsize=8, ncol=2)

    fig.set_size_inches(2, 1.5)
    fig.set_dpi(300)

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)

    plt.savefig(os.path.join(OUTPUT_FOLDER, "FigS2A.png"), bbox_inches='tight', dpi=300)

    # Fig S2B - scatter plot - Combfold vs. AFMv3
    MAX_T = 1
    combfold_tm_results, afm2_tm_results = [], []
    for jobname in pdb_to_subunits:
        combfold_tm_results.append(max(combfold_parsed_results.get(jobname, [])[:MAX_T], default=0))
        afm2_tm_results.append(max(afm3_results.get(jobname, [])[:MAX_T], default=0))
        if afm2_tm_results[-1] > 0.7 > combfold_tm_results[-1]:
            print(jobname, "better on AFMv3", afm2_tm_results[-1], combfold_tm_results[-1])

    fig, ax = plt.subplots(figsize=(2.25, 1.75), dpi=600)
    ax.scatter(combfold_tm_results, afm2_tm_results, alpha=0.4, s=10, edgecolor='none')
    ax.set_xlabel("CombFold", fontsize=11)
    ax.set_ylabel("AFMv3", fontsize=11)
    plt.yticks(np.arange(0, 1.2, 0.2), fontsize=8)
    plt.xticks(np.arange(0, 1.2, 0.2), fontsize=8)

    ax.text(0.5, 1.1, "Top-1 TM-score", ha='center', va='center', transform=ax.transAxes, fontsize=11)

    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.plot([0, 1], [0, 1], 'r--')

    plt.savefig(os.path.join(OUTPUT_FOLDER, "FigS2B.png"), bbox_inches='tight', dpi=300)


if __name__ == "__main__":
    main()