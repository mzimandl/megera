program=megera
OBJ=main.o cObjects.o cEvolver.o cDataGrid.o math_phys.o preset.o sSettings.o GenerateObj.o AnalyzeObj.o drawing.o
OPT=-O2 `freetype-config --cflags` -I ~/local/include -L ~/local/lib -lpng -lpngwriter -lfreetype -lz -fopenmp

${program}: ${OBJ}
	g++ ${OBJ} -o ${program} ${OPT}
	
main.o: main.cpp
	g++ main.cpp -c ${OPT}

cObjects.o: cObjects.cpp
	g++ cObjects.cpp -c ${OPT}

cEvolver.o: cEvolver.cpp
	g++ cEvolver.cpp -c ${OPT}

cDataGrid.o: cDataGrid.cpp
	g++ cDataGrid.cpp -c ${OPT}

math_phys.o: math_phys.cpp
	g++ math_phys.cpp -c ${OPT}

preset.o: preset.cpp
	g++ preset.cpp -c ${OPT}

sSettings.o: sSettings.cpp
	g++ sSettings.cpp -c ${OPT}

GenerateObj.o: GenerateObj.cpp
	g++ GenerateObj.cpp -c ${OPT}

AnalyzeObj.o: AnalyzeObj.cpp
	g++ AnalyzeObj.cpp -c ${OPT}

drawing.o: drawing.cpp
	g++ drawing.cpp -c ${OPT}