CC=g++
CFLAGS=-std=c++11 -O3
CFLAGS += $(shell $(CC) -fopenmp -E - < /dev/null > /dev/null 2>&1 && echo "-fopenmp")
LDFLAGS += $(shell $(CC) -fopenmp -E - < /dev/null > /dev/null 2>&1 && echo "-fopenmp")

.PHONY : clean linearalifold
objects= bin/linearalifold

linearalifold: src/Linearalifold.cpp
	mkdir -p bin
	$(CC) src/Linearalifold.cpp src/Utils/energy_model.cpp $(CFLAGS) -Dlv -Dis_candidate_list -o bin/linearalifold $(LDFLAGS)

clean:
	-rm $(objects)
