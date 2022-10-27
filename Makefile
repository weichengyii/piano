SOURCE = *.cpp
HEADER = *.h

piano: $(SOURCE) $(HEADER)
	g++ -g -o piano -O3 -funroll-loops -fomit-frame-pointer -falign-loops=16 $(SOURCE) -lm -lsndfile
