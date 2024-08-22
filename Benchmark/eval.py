import simexpal
import os
import pandas as pd
import numpy as np
import functools as ft

def aggregateParametricIBFS(df):
	df['iterationsPerBreakpoint'] = df['iterations']/df['breakpoints']
	df['bottlenecksPerIteration'] = df['bottlenecks']/df['iterations']
	df['bottlenecksPerBreakpoint'] = df['bottlenecks']/df['breakpoints']
	df['adoptionsPerBottleneck'] = df['adoptions']/df['bottlenecks']
	df['adoptionsPerIteration'] = df['adoptions']/df['iterations']
	df['adoptionsPerBreakpoint'] = df['adoptions']/df['breakpoints']
	df['loopToInitRatio'] = (df['updateTime'] + df['reconnectTime'] + df['drainTime'])/df['initTime']

def aggregateChordScheme(df):
	df['avgVertices'] = df['totalVertices']/df['vertices'] - 2
	df['contractionRatio'] = df['contractionTime']/df['runtime']

def aggregateSimexpal():
	cfg = simexpal.config_for_dir()

	dict = {}
	for successful_run in cfg.collect_successful_results():
		file = successful_run.output_file_path(ext="stats.csv")
		with open(file) as f:
			key = successful_run.experiment.name
			if not key in dict:
				dict[key] = []
			dict[key].append(pd.read_csv(file))

	frames = {}
	for key in dict:
		df = pd.concat(dict[key], ignore_index=True)
		if key == "parametricIBFS":
			aggregateParametricIBFS(df)
		elif key == "chordScheme":
			aggregateChordScheme(df)
		frames[key] = df
		df.to_csv("results/" + key + ".csv", index=False)

	return frames

def selectColumns(df, selector):
	return df.rename(columns=selector)[selector.values()]

def createMainTable(frames):
	selectorPBFS = {'instance': 'instance', 'vertices': 'vertices', 'edges': 'edges', 'runtime': 'runtimePBFS',
					'breakpoints': 'breakpointsPBFS', 'adoptionsPerBreakpoint': 'adoptionsPerBreakpointPBFS',
					'loopToInitRatio': 'loopToInitRatioPBFS'}
	pbfs = selectColumns(frames['parametricIBFS'], selectorPBFS)

	chordIBFS = frames['chordScheme'].loc[lambda df: df['algorithm'] == 'chordScheme[IBFS]']
	selectorChordIBFS = {'instance': 'instance', 'runtime': 'runtimeDSIBFS', 'breakpoints': 'breakpointsDSIBFS',
						 'avgVertices': 'avgVerticesDSIBFS', 'contractionRatio': 'contractionRatioDSIBFS'}
	chordIBFS = selectColumns(chordIBFS, selectorChordIBFS)

	tbl = pbfs.merge(chordIBFS, on='instance')
	tbl['speedup'] = tbl['runtimeDSIBFS']/tbl['runtimePBFS']
	tbl.to_csv("results/mainTable.csv", index=False)

def createAllAlgorithmsTable(frames):
	selectorPBFS = {'instance': 'instance', 'vertices': 'vertices', 'edges': 'edges', 'breakpoints': 'breakpoints',
					'runtime': 'runtimePBFS'}
	pbfs = selectColumns(frames['parametricIBFS'], selectorPBFS)

	chordIBFS = frames['chordScheme'].loc[lambda df: df['algorithm'] == 'chordScheme[IBFS]']
	selectorChordIBFS = {'instance': 'instance', 'runtime': 'runtimeDSIBFS'}
	chordIBFS = selectColumns(chordIBFS, selectorChordIBFS)

	chordPRF = frames['chordScheme'].loc[lambda df: df['algorithm'] == 'chordScheme[PushRelabel]']
	selectorChordPRF = {'instance': 'instance', 'runtime': 'runtimeDSPRF'}
	chordPRF = selectColumns(chordPRF, selectorChordPRF)

	tbl = pbfs.merge(chordIBFS, on='instance').merge(chordPRF, on='instance')
	tbl['speedupPBFSIBFS'] = tbl['runtimeDSIBFS']/tbl['runtimePBFS']
	tbl['speedupPBFSPRF'] = tbl['runtimeDSPRF']/tbl['runtimePBFS']
	tbl['speedupIBFSPRF'] = tbl['runtimeDSPRF']/tbl['runtimeDSIBFS']
	tbl.to_csv("results/allAlgorithmsTable.csv", index=False)

def createLiverTable(frames):
	selectorPBFS = {'instance': 'instance', 'runtime': 'runtimePBFS', 'breakpoints': 'breakpointsPBFS',
					'bottlenecksPerBreakpoint': 'bottlenecksPerBreakpointPBFS',
					'adoptionsPerBottleneck': 'adoptionsPerBottleneckPBFS',
					'avgDistance': 'avgDistancePBFS', 'loopToInitRatio': 'loopToInitRatioPBFS'}
	pbfs = selectColumns(frames['parametricIBFS'], selectorPBFS)

	chordIBFS = frames['chordScheme'].loc[lambda df: df['algorithm'] == 'chordScheme[IBFS]']
	selectorChordIBFS = {'instance': 'instance', 'runtime': 'runtimeDSIBFS', 'breakpoints': 'breakpointsDSIBFS',
						 'avgVertices': 'avgVerticesDSIBFS', 'contractionRatio': 'contractionRatioDSIBFS'}
	chordIBFS = selectColumns(chordIBFS, selectorChordIBFS)

	tbl = pbfs.merge(chordIBFS, on='instance')
	tbl = tbl.loc[tbl['instance'].str.contains('liver-n6c10-')]
	tbl['speedup'] = tbl['runtimeDSIBFS']/tbl['runtimePBFS']
	tbl.to_csv("results/liverTable.csv", index=False)

#TODO: Columns per instance
def createEpsilonTable(frames):
	selectorPBFS = {'instance': 'instance', 'runtime': 'runtimePBFS'}
	pbfs = selectColumns(frames['parametricIBFS'], selectorPBFS)

	selectorChord = {'instance': 'instance', 'epsilon': 'exponent', 'runtime': 'runtimeDS',
					 'breakpoints': 'breakpoints'}
	chord = selectColumns(frames['chordPrecision'], selectorChord)

	tbl = pbfs.merge(chord, on='instance')
	tbl['slowdown'] = tbl['runtimeDS']/tbl['runtimePBFS']
	tbl = tbl.drop(columns=['runtimePBFS', 'runtimeDS'])
	np.seterr(divide = 'ignore')
	tbl['exponent'] = np.abs(np.log10(tbl['exponent'])).replace(np.inf, np.nan).astype("Int32")
	np.seterr(divide = 'warn')
	tbl = tbl.sort_values(by=['exponent'])

	tbls = []
	for instance in tbl['instance'].unique():
		filtered = tbl.loc[tbl['instance'] == instance].copy()
		filtered['breakpoints'] = filtered['breakpoints']/filtered.iloc[-1]['breakpoints']
		selector = {"exponent": "exponent", "breakpoints": "breakpoints_" + instance, "slowdown": "slowdown_" + instance}
		filtered = selectColumns(filtered, selector)
		tbls.append(filtered)
	tbl = ft.reduce(lambda left, right: pd.merge(left=left, right=right, how='left', left_on=['exponent'],right_on=['exponent']), tbls)
	tbl.to_csv("results/epsilonTable.csv", index=False, sep='\t')

#TODO: Proper formatting of table entries
if __name__ == "__main__":
	if not os.path.exists("results/"):
		os.makedirs("results/")

	frames = aggregateSimexpal()
	createMainTable(frames)
	createAllAlgorithmsTable(frames)
	createLiverTable(frames)
	createEpsilonTable(frames)