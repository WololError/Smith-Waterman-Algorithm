all: projet

projetprelim: projetprelim.cpp
	g++ main.cpp src/* -o main

clean:
	rm *.o
