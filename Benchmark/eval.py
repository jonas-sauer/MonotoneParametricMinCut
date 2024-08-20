import simexpal
import os
import pandas as pd

if __name__ == "__main__":
	if not os.path.exists("results/"):
		os.makedirs("results/")

	cfg = simexpal.config_for_dir()
	dict = {}

	for successful_run in cfg.collect_successful_results():
		file = successful_run.output_file_path(ext="stats.csv")
		with open(file) as f:
			if not successful_run.experiment.name in dict:
				dict[successful_run.experiment.name] = []
			dict[successful_run.experiment.name].append(pd.read_csv(file))

	for key in dict:
		pd.concat(dict[key], ignore_index=True).to_csv("results/" + key + ".csv", index=False)
