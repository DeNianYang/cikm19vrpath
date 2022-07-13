OBJS = Tools.o LocoState.o SearchState.o SearchSpace.o Walk.o VirIndex.o PhyIndex.o DPCvsp.o kshortest.o mcpp.o
CC = g++ -O3 -std=c++11
CK = g++ -I.. -std=c++11
DEBUG = -g
CFLAGS = -c $(DEBUG)
LFLAGS = $(DEBUG)
all: CVSPsearch

CVSPsearch : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) main.cpp -o CVSPsearch

Tools.o : Tools.h Tools.cpp
	$(CC) $(CFLAGS) Tools.cpp

LocoState.o : LocoState.h LocoState.cpp Tools.h
	$(CC) $(CFLAGS) LocoState.cpp

Walk.o : Walk.h Walk.cpp Tools.h LocoState.h
	$(CC) $(CFLAGS) Walk.cpp

SearchState.o : SearchState.h SearchState.cpp LocoState.h Tools.h
	$(CC) $(CFLAGS) SearchState.cpp

SearchSpace.o : SearchSpace.h SearchSpace.cpp SearchState.h VirIndex.h PhyIndex.h
	$(CC) $(CFLAGS) SearchSpace.cpp

VirIndex.o : VirIndex.h VirIndex.cpp Tools.h
	$(CC) $(CFLAGS) VirIndex.cpp

PhyIndex.o : PhyIndex.h PhyIndex.cpp Tools.h Walk.h LocoState.h
	$(CC) $(CFLAGS) PhyIndex.cpp

DPCvsp.o : DPCvsp.h DPCvsp.cpp Tools.h Walk.h LocoState.h
	$(CC) $(CFLAGS) DPCvsp.cpp

kshortest.o: kshortestpath.cpp
	$(CK) $(CFLAGS) -c kshortestpath.cpp -o $@

mcpp.o: MCPP.cpp
	$(CK) $(CFLAGS) -c MCPP.cpp -o $@

clean:
	\rm *.o *~ a




