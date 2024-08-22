import simexpal
import os
import pandas as pd
import numpy as np
import functools as ft

def flatten(pivot):
	pivot.columns = pivot.columns.to_series().str.join("_")

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
	if not 'parametricIBFS' in frames:
		return
	if not 'chordScheme' in frames:
		return

	selectorPBFS = {'instance': 'instance', 'vertices': 'vertices', 'edges': 'edges', 'runtime': 'runtime_PBFS',
					'breakpoints': 'breakpoints_PBFS', 'adoptionsPerBreakpoint': 'adoptionsPerBreakpoint_PBFS',
					'loopToInitRatio': 'loopToInitRatio_PBFS'}
	pbfs = selectColumns(frames['parametricIBFS'], selectorPBFS)

	chordIBFS = frames['chordScheme'].loc[lambda df: df['algorithm'] == 'chordScheme[IBFS]']
	selectorChordIBFS = {'instance': 'instance', 'runtime': 'runtime_chordScheme[IBFS]', 'breakpoints': 'breakpoints_chordScheme[IBFS]',
						 'avgVertices': 'avgVertices_chordScheme[IBFS]', 'contractionRatio': 'contractionRatio_chordScheme[IBFS]'}
	chordIBFS = selectColumns(chordIBFS, selectorChordIBFS)

	tbl = pbfs.merge(chordIBFS, on='instance')
	tbl['speedup'] = (tbl['runtime_chordScheme[IBFS]']/tbl['runtime_PBFS']).round(2)
	tbl['runtime_PBFS'] = (tbl['runtime_PBFS']/1000).round(1)
	tbl['adoptionsPerBreakpoint_PBFS'] = tbl['adoptionsPerBreakpoint_PBFS'].round(1)
	tbl['loopToInitRatio_PBFS'] = tbl['loopToInitRatio_PBFS'].round(2)
	tbl['runtime_chordScheme[IBFS]'] = (tbl['runtime_chordScheme[IBFS]']/1000).round(1)
	tbl['avgVertices_chordScheme[IBFS]'] = tbl['avgVertices_chordScheme[IBFS]'].round(2)
	tbl['contractionRatio_chordScheme[IBFS]'] = (tbl['contractionRatio_chordScheme[IBFS]']*100).astype(int).astype(str) + "%"
	tbl.to_csv("results/mainTable.csv", index=False)

def createAllAlgorithmsTable(frames):
	if not 'parametricIBFS' in frames:
		return
	if not 'chordScheme' in frames:
		return

	selectorPBFS = {'instance': 'instance', 'vertices': 'vertices', 'edges': 'edges', 'breakpoints': 'breakpoints',
					'runtime': 'runtime_PBFS'}
	pbfs = selectColumns(frames['parametricIBFS'], selectorPBFS)

	chord = frames['chordScheme'].copy()
	chord['runtime'] = (chord['runtime']/1000).round(1)
	pivot = chord.pivot_table(index=['instance'], columns=['algorithm'], values=['runtime'])
	flatten(pivot)

	tbl = pbfs.merge(pivot, on='instance')
	tbl['runtime_PBFS'] = (tbl['runtime_PBFS']/1000).round(1)
	tbl['speedup_PBFS_IBFS'] = (tbl['runtime_chordScheme[IBFS]']/tbl['runtime_PBFS']).round(2)
	tbl['speedup_PBFS_PRF'] = (tbl['runtime_chordScheme[PushRelabel]']/tbl['runtime_PBFS']).round(2)
	tbl['speedup_IBFS_PRF'] = (tbl['runtime_chordScheme[PushRelabel]']/tbl['runtime_chordScheme[IBFS]']).round(2)
	tbl.to_csv("results/allAlgorithmsTable.csv", index=False)

def createLiverTable(frames):
	if not 'parametricIBFS' in frames:
		return
	if not 'chordScheme' in frames:
		return

	selectorPBFS = {'instance': 'instance', 'runtime': 'runtime_PBFS', 'breakpoints': 'breakpoints_PBFS',
					'bottlenecksPerBreakpoint': 'bottlenecksPerBreakpoint_PBFS',
					'adoptionsPerBottleneck': 'adoptionsPerBottleneck_PBFS',
					'avgDistance': 'avgDistance_PBFS', 'loopToInitRatio': 'loopToInitRatio_PBFS'}
	pbfs = selectColumns(frames['parametricIBFS'], selectorPBFS)

	chordIBFS = frames['chordScheme'].loc[lambda df: df['algorithm'] == 'chordScheme[IBFS]']
	selectorChordIBFS = {'instance': 'instance', 'runtime': 'runtime_chordScheme[IBFS]',
						 'breakpoints': 'breakpoints_chordScheme[IBFS]', 'avgVertices': 'avgVertices_chordScheme[IBFS]',
						 'contractionRatio': 'contractionRatio_chordScheme[IBFS]'}
	chordIBFS = selectColumns(chordIBFS, selectorChordIBFS)

	tbl = pbfs.merge(chordIBFS, on='instance')
	prefix = 'liver-n6c10-par-'
	tbl = tbl.loc[tbl['instance'].str.contains(prefix)]
	tbl['instance'] = tbl['instance'].str.replace(prefix, '')
	temp = tbl['instance'].str.split('-').str
	tbl.insert(loc=0, column='p_src', value=temp[0].astype(float)/10)
	tbl.insert(loc=1, column='p_snk', value=temp[1].astype(float)/10)
	tbl = tbl.drop('instance', axis=1)
	tbl['speedup'] = (tbl['runtime_chordScheme[IBFS]']/tbl['runtime_PBFS']).round(2)
	tbl['runtime_PBFS'] = (tbl['runtime_PBFS']/1000).round(1)
	tbl['bottlenecksPerBreakpoint_PBFS'] = tbl['bottlenecksPerBreakpoint_PBFS'].round(1)
	tbl['adoptionsPerBottleneck_PBFS'] = tbl['adoptionsPerBottleneck_PBFS'].round(1)
	tbl['avgDistance_PBFS'] = tbl['avgDistance_PBFS'].round(1)
	tbl['loopToInitRatio_PBFS'] = tbl['loopToInitRatio_PBFS'].round(2)
	tbl['runtime_chordScheme[IBFS]'] = (tbl['runtime_chordScheme[IBFS]']/1000).round(1)
	tbl['avgVertices_chordScheme[IBFS]'] = tbl['avgVertices_chordScheme[IBFS]'].round(2)
	tbl['contractionRatio_chordScheme[IBFS]'] = (tbl['contractionRatio_chordScheme[IBFS]']*100).astype(int).astype(str) + "%"
	tbl.to_csv("results/liverTable.csv", index=False)

def createEpsilonTable(frames):
	if not 'parametricIBFS' in frames:
		return
	if not 'chordPrecision' in frames:
		return

	selectorPBFS = {'instance': 'instance', 'runtime': 'runtime_PBFS'}
	pbfs = selectColumns(frames['parametricIBFS'], selectorPBFS)

	selectorChord = {'instance': 'instance', 'epsilon': 'exponent', 'runtime': 'runtime_DS',
					 'breakpoints': 'breakpoints'}
	chord = selectColumns(frames['chordPrecision'], selectorChord)

	tbl = pbfs.merge(chord, on='instance')
	tbl['slowdown'] = tbl['runtime_DS']/tbl['runtime_PBFS']
	tbl = tbl.drop(columns=['runtime_PBFS', 'runtime_DS'])
	np.seterr(divide = 'ignore')
	tbl['exponent'] = np.abs(np.log10(tbl['exponent'])).replace(np.inf, 1000).astype(int)
	np.seterr(divide = 'warn')
	tbl = tbl.sort_values(by=['exponent'])

	pivot = tbl.pivot_table(index=['exponent'], columns=['instance'], values=['breakpoints', 'slowdown'])
	pivot['breakpoints'] = pivot['breakpoints']/pivot['breakpoints'].iloc[-1]
	flatten(pivot)
	pivot.to_csv("results/epsilonTable.csv", index=False, sep='\t')

if __name__ == "__main__":
	if not os.path.exists("results/"):
		os.makedirs("results/")

	frames = aggregateSimexpal()
	createMainTable(frames)
	createAllAlgorithmsTable(frames)
	createLiverTable(frames)
	createEpsilonTable(frames)