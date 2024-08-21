import os

rawInstancesPath = "../Data/RawInstances"
instancesPath = "../Data/Instances"

def exec(cmd):
	print(cmd)
	os.system(cmd)

if __name__ == "__main__":
	for root, _, files in os.walk(rawInstancesPath):
		for file in files:
			if not file.endswith(".tbz2"):
				continue
			exec("tar -xvjf " + root + "/" + file + " -C " + root)