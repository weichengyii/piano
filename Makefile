SOURCE = piano.cpp hammer.cpp filter.cpp dwgs.cpp reverb.cpp sndfile.cpp
HEADER = types.h filter.h dwgs.h hammer.h reverb.h sndfile.h

piano: $(SOURCE) $(HEADER)
	g++ -g -o piano -O3 -funroll-loops -fomit-frame-pointer -falign-loops=16 $(SOURCE) -lm -lsndfile
