CC = g++
CFLAGS = -O3 -std=c++0x -fopenmp -funroll-loops -g -Wall
PRG = cycle_finder
OBJ = main.o common.o extract.o peak_detect.o kmer_compare.o cycle.o cycle_find.o trf_filter.o map_read.o blastn2fa.o repeat_num.o bll.o cluster.o kmer_align.o blast_to_copynum.o clustering.o clst_to_family.o intersperse.o path_find.o
DUMMY = common_dummy.h
DUMMY_SED = common.h



all: $(DUMMY_SED) $(PRG)

$(DUMMY_SED):
	@sed -e "s|ROOT_PATH_DUMMY|`pwd`|g" $(DUMMY) > $(DUMMY_SED)

$(PRG): $(OBJ)
	$(CC) -o $@ $(OBJ) $(CFLAGS)

.cpp.o:
	$(CC) -o $@ -c $< $(CFLAGS)

clean:
	rm -f $(PRG) $(OBJ) $(DUMMY_SED)