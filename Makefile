CFLAGS = -Wall -pedantic -O4
CC = g++
LDFLAGS = -Wl,-rpath
OBJECTS = functions.o ./CEC2013/cec2013.o ./CEC2013/cfunction.o ./CEC2013/rand2.o \
./Archive/Archive.o \
./BSP/RecordRevisitingScheme_Kernel/RecordRevisitingScheme_Kernel.o \
./BSP/VisitStrategy_Kernel/VisitStrategy.o \
./BSP/ACC_Math_Kernel/ACC_Math_Kernel.o \
./BSP/ClusterTree_Kernel/ClusterTree_Kernel.o \
./BSP/GeneticAlgorithm_Kernel/GeneticAlgorithm_Kernel.o \
./BSP/NonRevisitingScheme_Kernel/NonRevisitingScheme_Kernel.o \
./BSP/NonRevisitingGA_Kernel/NonRevisitingGA_Kernel.o \
./BSP/Optimization_TestFunction_Kernel/Optimization_TestFunction_Kernel.o \
./BSP/SelectionStrategey_Kernel/Selection_Strategy.o
INCLUDE = ./CEC2013/cec2013.h ./CEC2013/cfunction.h ./CEC2013/rand2.h  header.h \
./Archive/Archive.h \
./BSP/RecordRevisitingScheme_Kernel/RecordRevisitingScheme_Kernel.h \
./BSP/VisitStrategy_Kernel/VisitStrategy.h \
./BSP/ACC_Math_Kernel/ACC_Math_Kernel.h \
./BSP/ClusterTree_Kernel/ClusterTree_Kernel.h \
./BSP/GeneticAlgorithm_Kernel/GeneticAlgorithm_Kernel.h \
./BSP/NonRevisitingScheme_Kernel/NonRevisitingScheme_Kernel.h \
./BSP/NonRevisitingGA_Kernel/NonRevisitingGA_Kernel.h \
./BSP/Optimization_TestFunction_Kernel/Optimization_TestFunction_Kernel.h \
./BSP/SelectionStrategey_Kernel/Selection_Strategy.h


all: PNF_DE PNF_PSO

PNF_DE: PNF_DE.o $(OBJECTS) $(INCLUDE)
	$(CC) $(CFLAGS) -o PNF_DE PNF_DE.o $(OBJECTS)

PNF_PSO: PNF_PSO.o $(OBJECTS) $(INCLUDE)
	$(CC) $(CFLAGS) -o PNF_PSO PNF_PSO.o $(OBJECTS)

PNF_PSO.o: ./BSP/ACC_Math_Kernel/ACC_Math_Kernel.h \
	./BSP/ClusterTree_Kernel/ClusterTree_Kernel.h \
	./BSP/GeneticAlgorithm_Kernel/GeneticAlgorithm_Kernel.h \
	./BSP/NonRevisitingGA_Kernel/NonRevisitingGA_Kernel.h \
	./BSP/NonRevisitingScheme_Kernel/NonRevisitingScheme_Kernel.h \
	./BSP/Optimization_TestFunction_Kernel/Optimization_TestFunction_Kernel.h \
	./BSP/RecordRevisitingScheme_Kernel/RecordRevisitingScheme_Kernel.h \
	./BSP/VisitStrategy_Kernel/VisitStrategy.h \
	./BSP/SelectionStrategey_Kernel/Selection_Strategy.h 

PNF_DE.o: ./BSP/ACC_Math_Kernel/ACC_Math_Kernel.h \
	./BSP/ClusterTree_Kernel/ClusterTree_Kernel.h \
	./BSP/GeneticAlgorithm_Kernel/GeneticAlgorithm_Kernel.h \
	./BSP/NonRevisitingGA_Kernel/NonRevisitingGA_Kernel.h \
	./BSP/NonRevisitingScheme_Kernel/NonRevisitingScheme_Kernel.h \
	./BSP/Optimization_TestFunction_Kernel/Optimization_TestFunction_Kernel.h \
	./BSP/RecordRevisitingScheme_Kernel/RecordRevisitingScheme_Kernel.h \
	./BSP/VisitStrategy_Kernel/VisitStrategy.h \
	./BSP/SelectionStrategey_Kernel/Selection_Strategy.h 

functions.o: ./header.h ./BSP/ACC_Math_Kernel/ACC_Math_Kernel.h \
	./BSP/ClusterTree_Kernel/ClusterTree_Kernel.h \
	./BSP/GeneticAlgorithm_Kernel/GeneticAlgorithm_Kernel.h \
	./BSP/NonRevisitingGA_Kernel/NonRevisitingGA_Kernel.h \
	./BSP/NonRevisitingScheme_Kernel/NonRevisitingScheme_Kernel.h \
	./BSP/Optimization_TestFunction_Kernel/Optimization_TestFunction_Kernel.h \
	./BSP/VisitStrategy_Kernel/VisitStrategy.h


./BSP/VisitStrategy_Kernel/VisitStrategy.o: ./header.h \
	./Archive/Archive.h \
	./BSP/ACC_Math_Kernel/ACC_Math_Kernel.h \
	./BSP/GeneticAlgorithm_Kernel/GeneticAlgorithm_Kernel.h \
	./BSP/NonRevisitingGA_Kernel/NonRevisitingGA_Kernel.h \
	./BSP/Optimization_TestFunction_Kernel/Optimization_TestFunction_Kernel.h \
	./BSP/RecordRevisitingScheme_Kernel/RecordRevisitingScheme_Kernel.h \
	./BSP/SelectionStrategey_Kernel/Selection_Strategy.h 

./Archive/Archive.o: ./Archive/Archive.h

./BSP/RecordRevisitingScheme_Kernel/RecordRevisitingScheme_Kernel.o: ./BSP/ACC_Math_Kernel/ACC_Math_Kernel.h \
	./BSP/ClusterTree_Kernel/ClusterTree_Kernel.h \
	./BSP/GeneticAlgorithm_Kernel/GeneticAlgorithm_Kernel.h \
	./BSP/NonRevisitingScheme_Kernel/NonRevisitingScheme_Kernel.h \
	./BSP/SelectionStrategey_Kernel/Selection_Strategy.h 

./BSP/NonRevisitingGA_Kernel/NonRevisitingGA_Kernel.o: ./BSP/ACC_Math_Kernel/ACC_Math_Kernel.h \
	./BSP/ClusterTree_Kernel/ClusterTree_Kernel.h \
	./BSP/GeneticAlgorithm_Kernel/GeneticAlgorithm_Kernel.h \
	./BSP/NonRevisitingScheme_Kernel/NonRevisitingScheme_Kernel.h \
	./BSP/NonRevisitingGA_Kernel/NonRevisitingGA_Kernel.h

./BSP/NonRevisitingScheme_Kernel/NonRevisitingScheme_Kernel.o: ./BSP/NonRevisitingScheme_Kernel/NonRevisitingScheme_Kernel.h \
	./BSP/ClusterTree_Kernel/ClusterTree_Kernel.h

./BSP/GeneticAlgorithm_Kernel/GeneticAlgorithm_Kernel.o: ./BSP/ACC_Math_Kernel/ACC_Math_Kernel.h \
	./BSP/GeneticAlgorithm_Kernel/GeneticAlgorithm_Kernel.h

./BSP/SelectionStrategey_Kernel/Selection_Strategy.o: ./BSP/SelectionStrategey_Kernel/Selection_Strategy.h

./BSP/ClusterTree_Kernel/ClusterTree_Kernel.o:./BSP/ClusterTree_Kernel/ClusterTree_Kernel.h

./BSP/Optimization_TestFunction_Kernel/Optimization_TestFunction_Kernel.o: ./BSP/Optimization_TestFunction_Kernel/Optimization_TestFunction_Kernel.h \
	./BSP/ACC_Math_Kernel/ACC_Math_Kernel.h

./BSP/ACC_Math_Kernel/ACC_Math_Kernel.o: ./BSP/ACC_Math_Kernel/ACC_Math_Kernel.h
./CEC2013/cec2013.o: ./CEC2013/cfunction.h
./CEC2013/cfunction.o: ./CEC2013/cfunction.h
./CEC2013/rand2.o: ./CEC2013/rand2.h

.PhONY: clean
clean:
	rm -f $(OBJECTS) PNF_DE.o PNF_PSO.o