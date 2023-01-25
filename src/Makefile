myProgram: main.o initial_conds.o prim_and_con.o solvers.o outputFile.o constructor.o run.o header.H
	g++ -I/usr/include/eigen3 main.o initial_conds.o prim_and_con.o solvers.o outputFile.o constructor.o run.o -lconfig++ -o myProgram -O3

main.o: main.cpp
	g++ -c main.cpp -o main.o -O3

initial_conds.o: initial_conds.cpp
	g++ -c initial_conds.cpp -o initial_conds.o -O3

prim_and_con.o: prim_and_con.cpp
	g++ -c prim_and_con.cpp -o prim_and_con.o -O3

solvers.o: solvers.cpp
	g++ -c solvers.cpp -o solvers.o -O3

outputFile.o: outputFile.cpp
	g++ -c outputFile.cpp -o outputFile.o -O3

run.o: run.cpp
	g++ -c run.cpp -o run.o -O3

constructor.o: constructor.cpp
	g++ -c constructor.cpp -o constructor.o -O3
