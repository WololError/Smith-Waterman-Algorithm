all: projet

projetprelim: projetprelim.cpp
	g++ projetprelim.cpp src/* -o projetprelim

projet: projet.cpp
	g++ projet.cpp src/* -o projet

projetopt: projetopt.cpp
	g++ -O3 projetopt.cpp src/* -o projetopt

clean:
	rm *.o