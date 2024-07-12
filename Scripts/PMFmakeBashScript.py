import sys
from enum import Enum

configStuff = {"instanceSet": {}, "algorithmSet": {}, "modeSet": {}, "epsilonSet": {}}

runs = []

if __name__ == '__main__':
    configFileName: str
    if (sys.argv.__len__() <= 1):
        print("No config file given, reading PMFrunBenchmarkConfig.txt")
        configFileName = "PMFrunBenchmarkConfig.txt"
    else:
        configFileName = sys.argv[1]
        print("Reading config file " + configFileName)
    with open(configFileName, 'r') as config:
        currentSetName: str = ""
        currentSetType: str = ""
        runStuff = []
        runReadStep: int = 0
        for line in config.readlines():
            print(line)
            if str(line) == "\n" or line[0] == '#':
                continue
            if line[0] == '\\':
                words = line.split(' ')
                if str(line[1]).strip() == 'R':
                    print("Reading in run")
                    runStuff = []
                    runReadStep = 0
                    currentSetType = "run"
                    continue
                if str(line[1]).strip() == 'I':
                    print("Adding instance set" + words[1])
                    currentSetType = "instanceSet"
                if str(line[1]).strip() == 'A':
                    print("Adding algorithm set" + words[1])
                    currentSetType = "algorithmSet"
                if str(line[1]).strip() == 'M':
                    print("Adding mode set" + words[1])
                    currentSetType = "modeSet"
                if str(line[1]).strip() == 'E':
                    print("Adding epsilon set" + words[1])
                    currentSetType = "epsilonSet"
                currentSetName = str(words[1]).strip()
                configStuff[currentSetType][currentSetName] = set()
                continue
            if currentSetType != "run":
                print("Adding line " + line + " to set " + currentSetName + " of type", currentSetType)
                print(configStuff)
                configStuff[currentSetType][currentSetName].add(str(line).strip())
            else:
                runReadName: str = ""
                if runReadStep == 0:
                    runReadName = "instanceSet"
                if runReadStep == 1:
                    runReadName = "algorithmSet"
                if runReadStep == 2:
                    runReadName = "modeSet"
                if runReadStep == 3:
                    runReadName = "epsilonSet"
                if runReadStep == 4:
                    print("To many info for one run")

                runStuff.append(set())

                setsChosen = line.split(' ')
                for setName in setsChosen:
                    runStuff[runReadStep].update(configStuff[runReadName][str(setName).strip()])

                runReadStep += 1

                if runReadStep == 4:
                    runs.append((runStuff[0], runStuff[1], runStuff[2], runStuff[3]))

    print(configStuff)
    print(runs)

    runsSoFar = set()

    with open("PMFrunBenchmark.sh", 'w') as shellScript:
        shellScript.write("#!/bin/bash\n\n")
        shellScript.write("mkdir -p ../Data/Output\n\n")
        for run in runs:
            for instance in run[0]:
                for algorithm in run[1]:
                    for mode in run[2]:
                        for epsilon in run[3]:
                            if (instance, algorithm, mode, epsilon) in runsSoFar:
                                continue
                            if (mode == "whole"):
                                shellScript.write(
                                    "../build/BenchmarkParametricMaxFlow -i ../Data/" + instance + " -o ../Data/Output/parametricRuntimes.csv -a " + algorithm + " -m " + mode + " -e " + epsilon + "\n")
                            else:
                                shellScript.write(
                                    "../build/BenchmarkParametricMaxFlow -i ../Data/" + instance + " -o ../Data/Output/" + algorithm + "_" + mode + ".csv -a " + algorithm + " -m " + mode + " -e " + epsilon + "\n")
                            shellScript.write(
                                "echo \"finished run -i ../Data/" + instance + " -a " + algorithm + " -m " + mode + " -e " + epsilon + "\"\n")
                            runsSoFar.add((instance, algorithm, mode, epsilon))
