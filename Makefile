all: projet

projetprelim: projetprelim.cpp
	g++ projetprelim.cpp src/fasta.cpp src/blast.cpp src/Protein.cpp -o projetprelim
	./projetprelim query/P00533.fasta database/uniprot_sprot.fasta

projet: projet.cpp
	g++ projet.cpp src/fasta.cpp src/blast.cpp src/Protein.cpp src/SmithWaterman.cpp src/blosum.cpp -o projet

projetopt:
	@echo "Modifiez le fichier Makefile pour permettre la compilation de votre projet"

clean:
	rm *.o