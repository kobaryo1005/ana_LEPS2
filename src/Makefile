TARGET = ../bin/ana
ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)
CXXFLAGS = $(ROOTCFLAGS) -Wall
CXXLIBS = $(ROOTLIBS)
CC = g++

OBJS =  main.o Root.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) -g -m64 $(CXXLIBS) $(OBJS) -I../include -o $@
%.o : %.cc
	$(CC) $(CXXFLAGS) -c $< -I../include

clean:
	rm -f $(TARGET) $(OBJS)
