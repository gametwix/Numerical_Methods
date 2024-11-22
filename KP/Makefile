CC = g++
CFLAGS = -c
ADDLIB = -lm


all: libVolt.so

libVolt.o: volterra_solv.cpp
	$(CC) $(CFLAGS) -fpic volterra_solv.cpp -o libVolt.o 

libVolt.so: libVolt.o
	$(CC) -shared libVolt.o -o libVolt.so $(ADDLIB)
