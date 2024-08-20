import simexpal
import os

def write(file, content):
	f = open(file, "w")
	f.write(content)
	f.close()

if __name__ == "__main__":
	if not os.path.exists("results/"):
		os.makedirs("results/")

	cfg = simexpal.config_for_dir()
	resultPBFS = "algorithm,instance,vertices,edges,breakpoints,runtime,iterations,bottlenecks,adoptions,avgDistance,drains,initTime,updateTime,reconnectTime,drainTime\n"
	resultChord = "algorithm,instance,epsilon,vertices,edges,breakpoints,runtime,contractionTime,flowTime,totalVertices\n"
	resultChordNoContract = "algorithm,instance,epsilon,vertices,edges,epsilon,breakpoints,runtime\n"
	resultChordPrecision = "algorithm,instance,epsilon,vertices,edges,breakpoints,runtime,contractionTime,flowTime,totalVertices\n"
	resultRestartable = "algorithm,instance,epsilon,vertices,edges,breakpoints,runtime,updateTime,flowTime\n"

	for successful_run in cfg.collect_successful_results():
		with open(successful_run.output_file_path(ext="stats.csv")) as f:
			if successful_run.experiment.name == "parametricIBFS":
				resultPBFS += f.readline()
			elif successful_run.experiment.name == "chordScheme":
				resultChord += f.readline()
			elif successful_run.experiment.name == "chordSchemeNoContraction":
				resultChordNoContract += f.readline()
			elif successful_run.experiment.name == "chordPrecision":
				resultChordPrecision += f.readline()
			elif successful_run.experiment.name == "restartable":
				resultRestartable += f.readline()

	write("results/pbfs.csv", resultPBFS)
	write("results/chord.csv", resultChord)
	write("results/chord_nocontract.csv", resultChordNoContract)
	write("results/chord_precision.csv", resultChordPrecision)
	write("results/restartable.csv", resultRestartable)