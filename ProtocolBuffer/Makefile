CC=g++
#CC=clang
FLAGS=-std=c++17 -msse4 -fopenmp
OPTIMIZATION=-O3
DEBUG=-rdynamic -Werror -Wpedantic -pedantic-errors -Wall -Wextra -Wparentheses -Wno-sign-compare -Wfatal-errors
#DEBUG=-rdynamic -Werror -Wpedantic -pedantic-errors -Wall -Wextra -Wparentheses -Wno-sign-compare -Wno-missing-braces -Wno-literal-conversion -Wfatal-errors
RELEASE=-ftree-vectorize -Wfatal-errors -DNDEBUG

LIBRARY_PATH=-L /usr/local/lib
DYN_LIBRARIES=-lprotobuf -lcurl

all: ParseRealTimeGTFS
	
ParseRealTimeGTFS:
	$(CC) $(LIBRARY_PATH) $(FLAGS) $(OPTIMIZATION) $(DEBUG) -o ParseRealTimeGTFS ParseRealTimeGTFS.cpp Protocols/gtfs-realtime.pb.cc Protocols/gtfs-realtime-OVapi.pb.cc $(DYN_LIBRARIES) 
	
