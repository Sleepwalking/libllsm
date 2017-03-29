#CROSS = x86_64-w64-mingw32-
CC = $(CROSS)gcc
LINK = $(CROSS)gcc
AR = $(CROSS)ar
CFLAGS = -DFP_TYPE=float -Og -g -std=c99 -Wall -fPIC -fopenmp $(CFLAGSEXT)
ARFLAGS = -rv
OUT_DIR = ./build
OBJS = $(OUT_DIR)/math-funcs.o $(OUT_DIR)/llsm-level0.o $(OUT_DIR)/llsm-level1.o $(OUT_DIR)/envelope.o
LIBPYIN = external/libpyin
LIBGVPS = $(LIBPYIN)/external/libgvps

default: $(OUT_DIR)/libllsm.a
test: $(OUT_DIR)/llsm-test
	$(OUT_DIR)/llsm-test test/arctic_a0001.wav

$(OUT_DIR)/llsm-test: $(OUT_DIR)/libllsm.a $(OUT_DIR)/ciglet.o test/test.c \
    $(LIBGVPS)/build/libgvps.a $(LIBPYIN)/build/libpyin.a
	$(LINK) $(CFLAGS) -o $(OUT_DIR)/llsm-test \
	    test/test.c \
	    $(OUT_DIR)/libllsm.a $(OUT_DIR)/ciglet.o \
	    $(LIBPYIN)/build/libpyin.a \
	    $(LIBGVPS)/build/libgvps.a -lm -Wno-unused-result -fopenmp

$(LIBGVPS)/build/libgvps.a:
	cd $(LIBGVPS); mkdir -p build; make

$(LIBPYIN)/build/libpyin.a: $(LIBGVPS)//build/libgvps.a
	cd $(LIBPYIN); mkdir -p build; make

$(OUT_DIR)/libllsm.a: $(OBJS)
	$(AR) $(ARFLAGS) $(OUT_DIR)/libllsm.a $(OBJS) $(LIBS)
	@echo Done.

$(OUT_DIR)/math-funcs.o : math-funcs.c math-funcs.h
$(OUT_DIR)/llsm-level0.o : llsm-level0.c llsm.h envelope.h math-funcs.h
$(OUT_DIR)/llsm-level1.o : llsm-level1.c llsm.h envelope.h math-funcs.h
$(OUT_DIR)/envelope.o : envelope.c envelope.h math-funcs.h

$(OUT_DIR)/ciglet.o : external/ciglet/ciglet.c
	mkdir -p build
	$(CC) $(CFLAGS) -o $(OUT_DIR)/ciglet.o -c external/ciglet/ciglet.c

$(OUT_DIR)/%.o : %.c
	mkdir -p build
	$(CC) $(CFLAGS) -o $(OUT_DIR)/$*.o -c $*.c

install: $(OUT_DIR)/libllsm.a
	cp $(OUT_DIR)/libllsm.a /usr/lib/
	cp llsm.h /usr/include/
	@echo Done.

clean:
	@echo 'Removing all temporary binaries... '
	@rm -f $(OUT_DIR)/libllsm.a $(OUT_DIR)/*.o
	@echo Done.

clear:
	@echo 'Removing all temporary binaries... '
	@rm -f $(OUT_DIR)/libllsm.a $(OUT_DIR)/*.o
	@echo Done.
