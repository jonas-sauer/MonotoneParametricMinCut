import os

rawInstancesPath = "../Data/RawInstances"
instancesPath = "../Data/Instances"
instances = []

def findInstances():
	for root, _, files in os.walk(rawInstancesPath):
		newRoot = root.replace(rawInstancesPath, instancesPath)
		for file in files:
			if not file.endswith(".max"):
				continue
			oldPath = root + "/" + file
			newPath = newRoot + "/" + os.path.splitext(file)[0].replace('.', '-')
			instances.append((oldPath, newPath))

loader = "../cmake-build-release-wsl/InstanceLoader"
infinity = 999999

def exec(cmd):
	print(cmd)
	os.system(cmd)

def loadInstance(mode, input, output, p_src = "", p_snk = ""):
	cmd = loader + " -m " + mode + " -i " + input + " -o " + output + " -inf " + str(infinity)
	if p_src:
		cmd += " -p_src " + p_src
	if p_snk:
		cmd += " -p_snk " + p_snk
	exec(cmd)

def convertToBinary():
	print("Converting to binary format...")
	for (oldPath, newPath) in instances:
		loadInstance("static", oldPath, newPath)

default_source_probability = 1.0
all_source_probabilities = { 0.1, 0.5, 1.0 }
default_sink_probability = 0.0
all_sink_probabilities = { 0.0, 0.1, 0.5, 1.0 }
detailed_analysis_instances = { "liver.n6c10" }

def makeInstanceParametric(instance, p_src, p_snk):
	parametricInstance = instance + "-par-" + str(p_src).replace('.', '') + "-" + str(p_snk).replace('.', '')
	loadInstance("convert", instance, parametricInstance, str(p_src), str(p_snk))

def makeParametric():
	print("Making instances parametric...")
	for (_,instance) in instances:
		if any(i.replace('.', '-') in instance for i in detailed_analysis_instances):
			for p_src in all_source_probabilities:
				for p_snk in all_sink_probabilities:
					makeInstanceParametric(instance, p_src, p_snk)
		else:
			makeInstanceParametric(instance, default_source_probability, default_sink_probability)

if __name__ == "__main__":
	findInstances()
	convertToBinary()
	makeParametric()
