CC=g++
DEPS= src/Linearalifold.h src/Utils/energy_parameter.h src/Utils/feature_weight.h src/Utils/intl11.h src/Utils/intl21.h src/Utils/intl22.h src/Utils/utility_v.h src/Utils/utility.h src/Utils/array_utils.h
CFLAGS=-std=c++17 -O3
.PHONY : clean linearalifold
objects= bin/linearalifold bin/linearalifold_new

all: linearalifold

linearalifold: src/Linearalifold.cpp src/Linearalifold_new.cpp src/Utils/array_utils.cpp $(DEPS)
	mkdir -p bin
	$(CC) src/Linearalifold.cpp src/Utils/array_utils.cpp $(CFLAGS) -Dlv -Dis_candidate_list -o bin/linearalifold
	$(CC) src/Linearalifold_new.cpp src/Utils/array_utils.cpp $(CFLAGS) -Dlv -Dis_candidate_list -o bin/linearalifold_new

clean:
	-rm $(objects)
