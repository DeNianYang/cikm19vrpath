CC=g++
CXXFLAGS= -std=c++11
OBJS=cola.o
all: cola

cola: $(OBJS) 
	$(CC) $(CXXFLAGS) $(OBJS) -o cola

cola.o: cola.cpp
	$(CC) $(CXXFLAGS) -c cola.cpp -o $@

clean:
	rm -rf cola *.o

