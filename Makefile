paint: run
	gnuplot graph.gnu

run: compile
	./mke-v85 -s 100 -l

compile:
	g++ -std=c++11 my1.cpp -o mke-v85