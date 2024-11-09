paint: run
	gnuplot graph.gnu

run: compile
	./mke-v85 -s 1 -c

compile:
	g++ -std=c++11 mke-v85.cpp -o mke-v85