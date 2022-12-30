
CC = gcc
CFLAGS = -g -O2 -Wall -Wno-unused-function

all : ambimux

HTSDIR = htslib
include $(HTSDIR)/htslib.mk
include $(HTSDIR)/htslib_static.mk
HTSLIB = $(HTSDIR)/libhts.a

CPPFLAGS = -I. -I$(HTSDIR)

OBJS = main.o atac_data.o str_util.o sam_read.o \
	   variants.o bins.o overlap.o gtf_anno.o counts.o region.o \
	   clopts.o rna_data.o bc_stats.o mod.o bam_dat.o r_count.o \
	   bam_proc.o

LDFLAGS = -L$(HTSDIR)
LIBS = -lhts -lpthread $(HTSLIB_static_LIBS)

ambimux : $(OBJS) $(HTSLIB)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

%.o: %.c %.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $< 

clean:
	rm -f *o
	rm -f ambimux
	cd $(HTSDIR) && make clean && cd ../

