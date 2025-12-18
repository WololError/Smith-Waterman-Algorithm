all: projet

main: projetprelim.cpp
	g++ main.cpp src/* -o main

clean:
	rm *.o

