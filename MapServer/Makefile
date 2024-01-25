CC=g++
FLAGS=-O3 -std=c++17
DEBUG=-rdynamic -Werror -Wpedantic -pedantic-errors -Wall -Wextra -Wparentheses -Wno-sign-compare -Wfatal-errors

all: HalloWorldMapServer

HalloWorldMapServer:
	$(CC) $(FLAGS) -pthread -o mongoose.so -c mongoose/mongoose.c
	$(CC) $(FLAGS) $(DEBUG) -o HalloWorldMapServer HalloWorldMapServer.cpp -pthread mongoose.so -ldl
	rm mongoose.so
