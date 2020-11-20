CXX=g++
CXXFLAGS=-std=c++17 -Wall -g3
LDFLAGS=-lpqxx -lpq

AUX_OBJS=dbfunctions.o

MAIN_OBJS=$(AUX_OBJS) main.o

PROGS=main

all:	$(PROGS)

main:	$(MAIN_OBJS)
						$(CXX) $(CXXFLAGS) -o $@ $(MAIN_OBJS) $(LDFLAGS)
