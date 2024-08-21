import os

rawInstancesPath = "../Data/RawInstances"
instancesPath = "../Data/Instances"
loader = "../cmake-build-release-wsl/InstanceLoader"

def exec(cmd):
	print(cmd)
	os.system(cmd)

def loadInstance(input, output):
	cmd = loader + " -m parametric" + " -i " + input + " -o " + output
	exec(cmd)

if __name__ == "__main__":
	for root, _, files in os.walk(rawInstancesPath):
		newRoot = root.replace(rawInstancesPath, instancesPath)
		for file in files:
			if not file.endswith(".pmax"):
				continue
			oldPath = root + "/" + file
			newPath = newRoot + "/" + os.path.splitext(file)[0]
			loadInstance(oldPath, newPath)