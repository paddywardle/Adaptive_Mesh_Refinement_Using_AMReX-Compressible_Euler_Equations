EIGEN_PATH = "/usr/include/eigen3"

myProgram: main.o initial_conds.o prim_and_con.o solvers.o outputFile.o constructor.o run.o header.H
	g++ main.o initial_conds.o prim_and_con.o solvers.o outputFile.o constructor.o run.o -I /usr/include/eigen3 -o myProgram -O3

main.o: main.cpp header.H
	g++ -c main.cpp -o main.o -O3

initial_conds.o: initial_conds.cpp header.H
	g++ -c initial_conds.cpp -o initial_conds.o -O3

prim_and_con.o: prim_and_con.cpp header.H
	g++ -c prim_and_con.cpp -o prim_and_con.o -O3

solvers.o: solvers.cpp header.H
	g++ -c solvers.cpp -o solvers.o -O3

outputFile.o: outputFile.cpp header.H
	g++ -c outputFile.cpp -o outputFile.o -O3

run.o: run.cpp header.H
	g++ -c run.cpp -o run.o -O3

constructor.o: constructor.cpp header.H
	g++ -c constructor.cpp -o constructor.o -O3
