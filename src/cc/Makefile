all: main prof debug

main: main.cc
	g++ -O3 -std=c++11 -march=native main.cc -o pmc.opt
prof: main.cc
	g++ -g -pg -std=c++11 -march=native main.cc -o pmc.prf
debug: main.cc
	g++ -g -std=c++11 -march=native main.cc -o pmc.dbg


#	g++ -O3 -std=c++11 -march=native main.cc -o pmc
