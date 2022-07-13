CC=g++
CXXFLAGS= -std=c++11
OBJS=kernighanlin.o
all: kl

kl: $(OBJS) main.cpp
	$(CC) $(CXXFLAGS) $(OBJS) main.cpp -o kl

colapreprocessor.o: kernighanlin.cpp 
	$(CC) $(CXXFLAGS) -c kernighanlin.cpp -o $@

clean:
	rm -rf kl *.o

