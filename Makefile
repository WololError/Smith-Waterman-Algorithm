all: projet

projetprelim: projetprelim.cpp
	g++ projetprelim.cpp src/fasta.cpp src/blast.cpp src/Protein.cpp -o projetprelim

projet: projet.cpp
	g++ projet.cpp src/fasta.cpp src/blast.cpp src/Protein.cpp src/SmithWaterman.cpp src/blosum.cpp -o projet

projetopt: projetopt.cpp
	g++ -O3 -march=native -mtune=native -funroll-loops -ffast-math \
    projetopt.cpp src/fasta.cpp src/blast.cpp src/Protein.cpp \
    src/SmithWaterman.cpp src/blosum.cpp -o projet

clean:
	rm *.o