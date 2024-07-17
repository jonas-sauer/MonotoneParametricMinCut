NAME=CapacityAssignment
CC=g++ -fopenmp
FLAGS=-std=c++23 -pipe
OPTIMIZATION=-march=native -O3
RELEASE=-ffast-math -Wfatal-errors -DNDEBUG

INCLUDE=-I /usr/include/cairo
LIBRARY_PATH=-L /usr/lib
DYN_LIBRARIES=-lcairo

Build:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(RELEASE) -o Build Build.cpp $(DYN_LIBRARIES)

RunCommand:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(RELEASE) -o $(NAME).exe RunCommand.cpp $(DYN_LIBRARIES)
