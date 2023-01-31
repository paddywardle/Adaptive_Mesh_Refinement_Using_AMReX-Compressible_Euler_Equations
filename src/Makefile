myProgram: main.o EulerEOS.o EulerEOS.h EulerSystem.o EulerSystem.h NumericalMethod.o NumericalMethod.h SettingsData.o
	g++ -I/usr/include/eigen3 main.o EulerEOS.o EulerSystem.o NumericalMethod.o SettingsData.o -lconfig++ -o myProgram -O3

main.o: main.cpp
	g++ -c main.cpp -o main.o -O3

EulerEOS.o: EulerEOS.cpp EulerEOS.h
	g++ -c EulerEOS.cpp -o EulerEOS.o -O3

EulerSystem.o: EulerSystem.cpp EulerSystem.h
	g++ -c EulerSystem.cpp -o EulerSystem.o -O3

NumericalMethod.o: NumericalMethod.cpp NumericalMethod.h
	g++ -c NumericalMethod.cpp -o NumericalMethod.o -O3

SettingsData.o: SettingsData.cpp SettingsData.h
	g++ -c SettingsData.cpp -o SettingsData.o -O3
