NAME=CapacityAssignment
CC=g++ -fopenmp
FLAGS=-std=c++17 -pipe
OPTIMIZATION=-march=native -O3
DEBUG=-rdynamic -Werror -Wpedantic -pedantic-errors -Wall -Wextra -Wparentheses -Wfatal-errors -D_GLIBCXX_DEBUG -g -fno-omit-frame-pointer
RELEASE=-ffast-math -ftree-vectorize -Wfatal-errors -DNDEBUG

INCLUDE=-I /usr/include/cairo
LIBRARY_PATH=-L /usr/lib
DYN_LIBRARIES=-lcairo

Build:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(DEBUG) -o Build Build.cpp $(DYN_LIBRARIES)

RunCommand:
	$(CC) $(INCLUDE) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(RELEASE) -o $(NAME).exe RunCommand.cpp $(DYN_LIBRARIES)
