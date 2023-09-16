import json
import math
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator

THIS_SCRIPT_PATH = os.path.abspath(__file__)
DATA_PATH = os.path.join(os.path.dirname(THIS_SCRIPT_PATH), "data")
OUTPUT_FOLDER = os.path.join(os.path.dirname(THIS_SCRIPT_PATH), "output", "figureS7")


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

    def format_func(value, tick_number=None):
        num_thousands = 0 if abs(value) < 1000 else math.floor(math.log10(abs(value)) / 3)
        value = round(value / 1000 ** num_thousands, 2)
        return f'{value:g}' + ' KMGTPEZY'[num_thousands]
    ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(format_func))

    equation = f'$\\rho = {np.corrcoef(x, y)[0][1]:.2f}$'
    plt.annotate(equation, xy=(0.05, 0.9), xycoords='axes fraction', fontsize=8,
                 bbox=dict(boxstyle='square', facecolor='white', edgecolor="black", lw=0.5, pad=0.2))

    plt.savefig(output_path, bbox_inches='tight', dpi=300)


def get_duration_json(json_path):
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    json_data = json.load(open(json_path, "r"))
    without_dups = {}

    for jobname, times in json_data.items():
        if "multimer" not in jobname:
            continue
        if len(times) == 3 and times[0] == "dup":
            if len(times[2]) == 3 and times[2][0] == "dup":
                if len(times[2][2]) == 3 and times[2][2][0] == "dup":
                    without_dups[jobname] = times[1] + times[2][1] + times[2][2][1] + times[2][2][2]
                else:
                    without_dups[jobname] = times[1] + times[2][1] + times[2][2]
            else:
                without_dups[jobname] = times[1] + times[2]
        else:
            without_dups[jobname] = times

    parsed_times = {}
    for jobname, times in without_dups.items():
        if "dup" in times:
            print(jobname, times)
        parsed_times[jobname] = [float(i[:-1]) for i in times]
    return parsed_times


def main():

    # load benchmark info
    pdb_to_subunits = {}
    benchmark_path = os.path.join(DATA_PATH, "benchmark1")
    runtimes_folder = os.path.join(benchmark_path, "runtimes")

    for filename in os.listdir(os.path.join(benchmark_path, "pdb_infos")):
        jobname = filename.split(".")[0]
        if not filename.endswith(".json"):
            continue
        pdb_to_subunits[jobname] = json.load(open(os.path.join(benchmark_path, "pdb_infos", filename), "rb"))

    print("Benchmark size ", len(pdb_to_subunits), "PDB ids: ", list(pdb_to_subunits.keys()))

    combfold_results = json.load(open(os.path.join(benchmark_path, "combfold_results.json"), "r"))
    combfold_results = {k: v for k, v in combfold_results.items() if v is not None and k in pdb_to_subunits}
    combfold_assembly_time = {jobname: round(result["took"], 1) for jobname, result in combfold_results.items()}

    simple_afm_durations = get_duration_json(os.path.join(runtimes_folder, "duration_simple.json"))
    combfold_afm_durations = get_duration_json(os.path.join(runtimes_folder, "duration_combfold.json"))

    merged_afm_durations = {}
    for jobname, times in simple_afm_durations.items():
        jobname = jobname.split("_")[1]
        if jobname not in merged_afm_durations:
            merged_afm_durations[jobname] = []
        merged_afm_durations[jobname] += times
    avg_simple_afm_durations = {k: np.mean(v) for k, v in merged_afm_durations.items() if v}
    # print(avg_simple_afm_durations)
    print("Average AFM", np.mean(list(avg_simple_afm_durations.values())))

    avg_combfold_pairs_durations = {k: np.mean(v) for k, v in combfold_afm_durations.items()
                                    if v and len(k.split("_")) == 4}
    avg_combfold_groupss_durations = {k: np.mean(v) for k, v in combfold_afm_durations.items()
                                      if v and len(k.split("_")) != 4}

    print("Average CombFold pairs", np.mean(list(avg_combfold_pairs_durations.values())))
    print("Average CombFold groups", np.mean(list(avg_combfold_groupss_durations.values())))

    combfold_pairs_durations_by_jobname = {}
    for af_jobname, avg_model_time in avg_combfold_pairs_durations.items():
        jobname = af_jobname.split("_")[1]
        if jobname not in combfold_pairs_durations_by_jobname:
            combfold_pairs_durations_by_jobname[jobname] = []
        combfold_pairs_durations_by_jobname[jobname].append(avg_model_time)

    combfold_groups_durations_by_jobname = {}
    for af_jobname, avg_model_time in avg_combfold_groupss_durations.items():
        jobname = af_jobname.split("_")[1]
        if jobname not in combfold_groups_durations_by_jobname:
            combfold_groups_durations_by_jobname[jobname] = []
        combfold_groups_durations_by_jobname[jobname].append(avg_model_time)

    avg_jobname_combfold_pairs_durations = {k: np.mean(v) for k, v in combfold_pairs_durations_by_jobname.items()}
    avg_jobname_combfold_groups_durations = {k: np.mean(v) for k, v in combfold_groups_durations_by_jobname.items()}

    print("Average CombFold pairs by jobname", np.mean(list(avg_jobname_combfold_pairs_durations.values())))
    print("Average CombFold groups by jobname", np.mean(list(avg_jobname_combfold_groups_durations.values())))

    total_combfold_pairs_by_jobname = {}
    for af_jobname, avg_model_time in avg_combfold_pairs_durations.items():
        jobname = af_jobname.split("_")[1]
        if jobname not in total_combfold_pairs_by_jobname:
            total_combfold_pairs_by_jobname[jobname] = 0
        total_combfold_pairs_by_jobname[jobname] += avg_model_time
    print("Average Total CombFold pairs", np.mean(list(total_combfold_pairs_by_jobname.values())))

    total_combfold_groups_by_jobname = {}
    for af_jobname, avg_model_time in avg_combfold_groupss_durations.items():
        jobname = af_jobname.split("_")[1]
        if jobname not in total_combfold_groups_by_jobname:
            total_combfold_groups_by_jobname[jobname] = 0
        total_combfold_groups_by_jobname[jobname] += avg_model_time
    print("Average Total CombFold groups", np.mean(list(total_combfold_groups_by_jobname.values())))

    with open(os.path.join(runtimes_folder, "runtime_summary.csv"), "w") as f:
        f.write("PDB ID, #subunits, #unique subunits, AFM runtime, CombFold pairs average, CombFold groups average, "
                "CombFold pairs total, CombFold Groups total, CombFold assembly runtime\n")
        for jobname in sorted(list(pdb_to_subunits.keys())):
            total_chain_num = sum([len(i["chain_names"]) for i in pdb_to_subunits[jobname].values()])

            f.write(f"{jobname}, {total_chain_num}, {len(pdb_to_subunits[jobname])}, "
                    f"{avg_simple_afm_durations.get(jobname, '-')}, "
                    f"{avg_jobname_combfold_pairs_durations.get(jobname, '-')}, "
                    f"{avg_jobname_combfold_groups_durations.get(jobname, '-')}, "
                    f"{total_combfold_pairs_by_jobname.get(jobname, '-')}, "
                    f"{total_combfold_groups_by_jobname.get(jobname, '-')}, "
                    f"{combfold_assembly_time.get(jobname, '-')}\n")

    unique_subunits = []
    total_combfold_time = []
    for jobname in sorted(list(pdb_to_subunits.keys())):
        unique_subunits.append(len(pdb_to_subunits[jobname]))
        total_combfold_time.append(total_combfold_pairs_by_jobname.get(jobname, 0)
                                   + total_combfold_groups_by_jobname.get(jobname, 0)
                                   + combfold_assembly_time.get(jobname, 0))
    draw_scatter("Unique subunits", "Runtime (sec)", unique_subunits, total_combfold_time,
                 os.path.join(OUTPUT_FOLDER, "FigS7.png"))


if __name__ == "__main__":
    main()