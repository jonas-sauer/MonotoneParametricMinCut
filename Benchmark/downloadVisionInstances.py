import requests

url_path = "https://vision.cs.uwaterloo.ca/files/"
local_path = "../Data/RawInstances/Vision/"
prefixes = {"liver", "adhead", "babyface", "bone"}
configurations = {"n6c10", "n6c100", "n26c10", "n26c100"}

def download(file):
	print("Downloading " + file)
	resp = requests.get(url_path + file)
	with open(local_path + file, "wb") as f:
		f.write(resp.content)

if __name__ == "__main__":
	for prefix in prefixes:
		for config in configurations:
			download(prefix + '.' + config + '.tbz2')