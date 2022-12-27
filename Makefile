
CC = gcc
# CFLAGS = -g -O2 -Wall -Wno-unused-function -pg
CFLAGS = -g -O2 -Wall -Wno-unused-function

all : ambimux

HTSDIR = htslib
include $(HTSDIR)/htslib.mk
include $(HTSDIR)/htslib_static.mk
HTSLIB = $(HTSDIR)/libhts.a
HTSLIB_LIB = $(HTSLIB) $(HTSLIB_static_LIBS)
SCLIBS = $(HTSLIB_LIB) -lpthread

HTSLIB_CPPFLAGS = -I$(HTSDIR)
CPPFLAGS = -I. $(HTSLIB_CPPFLAGS)

OBJS = main.o atac_data.o str_util.o sam_read.o \
	   variants.o bins.o overlap.o gtf_anno.o counts.o region.o \
	   clopts.o rna_data.o bc_stats.o mod.o bam_dat.o r_count.o \
	   bam_proc.o

LDFLAGS = $(HTSLIB_LIBS)

ambimux : $(OBJS) $(HTSLIB)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(SCLIBS)

%.o: %.c %.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $< 

# for specific configuration but "include htslib.mk" seems to work well
# HTSCONF = --disable-lzma --without-curses
# $(HTSLIB) : 
# 	cd $(HTSDIR) && autoreconf -i && ./configure $(HTSCONF) && make && cd ../

clean:
	rm -f *o
	rm ambimux
	cd $(HTSDIR) && make clean && cd ../

