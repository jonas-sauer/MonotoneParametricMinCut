import requests
import os

filename = "PolygonAggregationInstances.zip"
url_path = "https://zenodo.org/records/13642985/files/" + filename + "?download=1"
local_path = "../Data/RawInstances/Aggregation/"
local_file = local_path + filename

if __name__ == "__main__":
	if not os.path.exists(local_path):
		os.makedirs(local_path)
	print("Downloading " + url_path)
	resp = requests.get(url_path)
	with open(local_file, "wb") as f:
		f.write(resp.content)
	os.system("unzip " + local_file + " -d " + local_path)