ideas:	main.cpp tensorHMMbase.cpp tensorHMMbase.h genomicTensor.cpp genomicTensor.h MixGauss.cpp MixGauss.h datastructure.h
	g++ -I. -I/gpfs/home/yzz2/work/gsl/include -c main.cpp tensorHMMbase.cpp genomicTensor.cpp MixGauss.cpp
	g++ -L/gpfs/home/yzz2/work/gsl/lib main.o tensorHMMbase.o genomicTensor.o MixGauss.o -lgsl -lgslcblas -lm -o ideas
