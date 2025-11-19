all: projet

projetprelim: projetprelim.cpp
	g++ projetprelim.cpp src/fasta.cpp src/blast.cpp src/Protein.cpp -o projetprelim
	./projetprelim query/P00533.fasta database/uniprot_sprot.fasta

projet: projet.cpp
	g++ projet.cpp src/*.cpp -o projet
	./projet query/P00533.fasta database/uniprot_sprot.fasta blosum/BLOSUM62 11 1

projetopt:
	@echo "Modifiez le fichier Makefile pour permettre la compilation de votre projet"

clean: 
	rm *.o
