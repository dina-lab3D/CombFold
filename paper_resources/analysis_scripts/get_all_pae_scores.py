import json
import os
import sys


def save_pae_from_benchmark(benchmark_path: str, output_folder: str):
    all_pae_scores = []
    used_pae_scores = []
    top1_used_pae_scores = []
    for jobname in os.listdir(benchmark_path):
        if not os.path.isdir(os.path.join(benchmark_path, jobname)):
            continue
        job_folder = os.path.join(benchmark_path, jobname)

        output_clusters_path = os.path.join(job_folder, "assembly_output", "output_clustered.res")
        if not os.path.exists(output_clusters_path):
            continue
        top1_score = 0
        top1_pae_scores = []
        for line in open(output_clusters_path, "r").readlines():
            if "foldSteps: " not in line:
                continue
            # print(line)
            # print(line.split("foldSteps: ")[1].split("-"))
            pae_scores = [float(i.split(" ")[0]) for i in line.split("foldSteps: ")[1].split("-")[1:] if i]
            used_pae_scores += pae_scores
            score = float(line.split(" ")[3])
            if score > top1_score:
                top1_score = score
                top1_pae_scores = pae_scores
        top1_used_pae_scores += top1_pae_scores

        for trans_file in os.listdir(os.path.join(job_folder, "transformations")):
            if not os.path.isfile(os.path.join(job_folder, "transformations", trans_file)):
                continue
            all_pae_scores += [float(i.split(" | ")[1]) for i in
                               open(os.path.join(job_folder, "transformations", trans_file), "r").readlines()
                               if " | " in i]
        print("Done jobname: ", jobname)

    json.dump(used_pae_scores, open(os.path.join(output_folder, "used_pae_scores.json"), "w"))
    json.dump(top1_used_pae_scores, open(os.path.join(output_folder, "top1_used_pae_scores.json"), "w"))
    json.dump(all_pae_scores, open(os.path.join(output_folder, "all_pae_scores.json"), "w"))


if __name__ == '__main__':
    assert len(sys.argv) == 3
    save_pae_from_benchmark(os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2]))
