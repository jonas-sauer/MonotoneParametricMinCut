CC=g++ -fopenmp
#CC=clang++-6.0
FLAGS=-std=c++17 -pipe
OPTIMIZATION=-march=native
DEBUG=-rdynamic -Werror -Wpedantic -pedantic-errors -Wall -Wextra -Wparentheses -Wfatal-errors -D_GLIBCXX_DEBUG -g -fno-omit-frame-pointer -O0
RELEASE=-ffast-math -ftree-vectorize -Wfatal-errors -DNDEBUG -O3

INCLUDE=-I /usr/include/cairo
LIBRARY_PATH=-L /usr/lib
DYN_LIBRARIES=-lcairo

allDebug: ShellsDebug MapServersDebug UnitTestsDebug

allRelease: ShellsRelease MapServersRelease UnitTestsRelease

#Shells
ShellsDebug: AssignmentDebug BenchmarkShellDebug BikeSharingDebug CHDebug DelayExperimentsShellDebug IsochroneShellDebug MultimodalCoverageShellDebug NetworkBuilderDebug NetworkOrderShellDebug PartitionShellDebug QueriesShellDebug ULTRADebug UltraAssignmentDebug ULTRAExperimentsShellDebug UnrestrictedWalkingShellDebug ValidationShellDebug

ShellsRelease: AssignmentRelease BenchmarkShellRelease BikeSharingRelease CHRelease DelayExperimentsShellRelease IsochroneShellRelease MultimodalCoverageShellRelease NetworkBuilderRelease NetworkOrderShellRelease PartitionShellRelease QueriesShellRelease ULTRARelease UltraAssignmentRelease ULTRAExperimentsShellRelease UnrestrictedWalkingShellRelease ValidationShellRelease

AssignmentDebug:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(DEBUG) -o Assignment Assignment.cpp $(DYN_LIBRARIES)

AssignmentRelease:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(RELEASE) -o Assignment Assignment.cpp $(DYN_LIBRARIES)

BenchmarkShellDebug:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(DEBUG) -o BenchmarkShell BenchmarkShell.cpp $(DYN_LIBRARIES)

BenchmarkShellRelease:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(RELEASE) -o BenchmarkShell BenchmarkShell.cpp $(DYN_LIBRARIES)

BikeSharingDebug:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(DEBUG) -o BikeSharing BikeSharing.cpp $(DYN_LIBRARIES)

BikeSharingRelease:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(RELEASE) -o BikeSharing BikeSharing.cpp $(DYN_LIBRARIES)

CHDebug:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(DEBUG) -o CH CH.cpp $(DYN_LIBRARIES)

CHRelease:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(RELEASE) -o CH CH.cpp $(DYN_LIBRARIES)

DelayExperimentsShellDebug:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(DEBUG) -o DelayExperimentsShell DelayExperimentsShell.cpp $(DYN_LIBRARIES)

DelayExperimentsShellRelease:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(RELEASE) -o DelayExperimentsShell DelayExperimentsShell.cpp $(DYN_LIBRARIES)

IsochroneShellDebug:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(DEBUG) -o IsochroneShell IsochroneShell.cpp $(DYN_LIBRARIES)

IsochroneShellRelease:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(RELEASE) -o IsochroneShell IsochroneShell.cpp $(DYN_LIBRARIES)

MultimodalCoverageShellDebug:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(DEBUG) -o MultimodalCoverageShell MultimodalCoverageShell.cpp $(DYN_LIBRARIES)

MultimodalCoverageShellRelease:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(RELEASE) -o MultimodalCoverageShell MultimodalCoverageShell.cpp $(DYN_LIBRARIES)

NetworkBuilderDebug:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(DEBUG) -o NetworkBuilder NetworkBuilder.cpp $(DYN_LIBRARIES)

NetworkBuilderRelease:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(RELEASE) -o NetworkBuilder NetworkBuilder.cpp $(DYN_LIBRARIES)

NetworkOrderShellDebug:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(DEBUG) -o NetworkOrderShell NetworkOrderShell.cpp $(DYN_LIBRARIES)

NetworkOrderShellRelease:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(RELEASE) -o NetworkOrderShell NetworkOrderShell.cpp $(DYN_LIBRARIES)

PartitionShellDebug:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(DEBUG) -o PartitionShell PartitionShell.cpp $(DYN_LIBRARIES)

PartitionShellRelease:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(RELEASE) -o PartitionShell PartitionShell.cpp $(DYN_LIBRARIES)

QueriesShellDebug:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(DEBUG) -o QueriesShell QueriesShell.cpp $(DYN_LIBRARIES)

QueriesShellRelease:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(RELEASE) -o QueriesShell QueriesShell.cpp $(DYN_LIBRARIES)

ULTRADebug:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(DEBUG) -o ULTRA ULTRA.cpp $(DYN_LIBRARIES)

ULTRARelease:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(RELEASE) -o ULTRA ULTRA.cpp $(DYN_LIBRARIES)

UltraAssignmentDebug:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(DEBUG) -o UltraAssignment UltraAssignment.cpp $(DYN_LIBRARIES)

UltraAssignmentRelease:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(RELEASE) -o UltraAssignment UltraAssignment.cpp $(DYN_LIBRARIES)

ULTRAExperimentsShellDebug:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(DEBUG) -o ULTRAExperimentsShell ULTRAExperimentsShell.cpp $(DYN_LIBRARIES)

ULTRAExperimentsShellRelease:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(RELEASE) -o ULTRAExperimentsShell ULTRAExperimentsShell.cpp $(DYN_LIBRARIES)

UnrestrictedWalkingShellDebug:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(DEBUG) -o UnrestrictedWalkingShell UnrestrictedWalkingShell.cpp $(DYN_LIBRARIES)

UnrestrictedWalkingShellRelease:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(RELEASE) -o UnrestrictedWalkingShell UnrestrictedWalkingShell.cpp $(DYN_LIBRARIES)

ValidationShellDebug:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(DEBUG) -o ValidationShell ValidationShell.cpp $(DYN_LIBRARIES)

ValidationShellRelease:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(RELEASE) -o ValidationShell ValidationShell.cpp $(DYN_LIBRARIES)

#MapServers
MapServersDebug: DemandMapServerDebug RaptorMapServerDebug

MapServersRelease: DemandMapServerRelease RaptorMapServerRelease

DemandMapServerDebug:
	$(CC) -std=c++14 -fopenmp -Wno-deprecated-register $(OPTIMIZATION) $(RELEASE) -pthread -o mongoose.so -c ../MapServer/mongoose/mongoose.c
	$(CC) $(FLAGS) $(INCLUDE) $(LIBRARY_PATH) $(OPTIMIZATION) $(DEBUG) -o DemandMapServer DemandMapServer.cpp -pthread mongoose.so -ldl $(DYN_LIBRARIES)
	rm mongoose.so

DemandMapServerRelease:
	$(CC) -std=c++14 -fopenmp -Wno-deprecated-register $(OPTIMIZATION) $(RELEASE) -pthread -o mongoose.so -c ../MapServer/mongoose/mongoose.c
	$(CC) $(FLAGS) $(INCLUDE) $(LIBRARY_PATH) $(OPTIMIZATION) $(RELEASE) -o DemandMapServer DemandMapServer.cpp -pthread mongoose.so -ldl $(DYN_LIBRARIES)
	rm mongoose.so

RaptorMapServerDebug:
	$(CC) -std=c++14 -fopenmp -Wno-deprecated-register $(OPTIMIZATION) $(RELEASE) -pthread -o mongoose.so -c ../MapServer/mongoose/mongoose.c
	$(CC) $(FLAGS) $(INCLUDE) $(LIBRARY_PATH) $(OPTIMIZATION) $(DEBUG) -o RaptorMapServer RaptorMapServer.cpp -pthread mongoose.so -ldl $(DYN_LIBRARIES)
	rm mongoose.so

RaptorMapServerRelease:
	$(CC) -std=c++14 -fopenmp -Wno-deprecated-register $(OPTIMIZATION) $(RELEASE) -pthread -o mongoose.so -c ../MapServer/mongoose/mongoose.c
	$(CC) $(FLAGS) $(INCLUDE) $(LIBRARY_PATH) $(OPTIMIZATION) $(RELEASE) -o RaptorMapServer RaptorMapServer.cpp -pthread mongoose.so -ldl $(DYN_LIBRARIES)
	rm mongoose.so

#Tests
UnitTestsDebug:
	$(CC) $(FLAGS) $(DEBUG) -o UnitTests UnitTests.cpp

UnitTestsRelease:
	$(CC) $(FLAGS) $(RELEASE) -o UnitTests UnitTests.cpp