
CC = gcc
# CFLAGS = -g -O1 -Wall -Wextra -Wfloat-equal -Wno-unused-function -fsanitize=address -Wshadow -fno-omit-frame-pointer
CFLAGS = -g -O2 -Wall -Wextra -Wfloat-equal -Wno-unused-function -Wpointer-arith -Wshadow

all : ambimux

HTSDIR = htslib
include $(HTSDIR)/htslib.mk
include $(HTSDIR)/htslib_static.mk
HTSLIB = $(HTSDIR)/libhts.a

CPPFLAGS = -I. -I$(HTSDIR)

OBJS = main.o atac_data.o str_util.o sam_read.o \
	   variants.o bins.o overlap.o gtf_anno.o counts.o region.o \
	   clopts.o rna_data.o bc_stats.o mod.o bam_dat.o r_count.o \
	   bam_proc.o math_util.o array_util.o

LDFLAGS = -L$(HTSDIR)
LIBS = -lm -lhts -lpthread $(HTSLIB_static_LIBS)

ambimux : $(HTSLIB) $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

%.o: %.c %.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $< 

cleano :
	rm -f *o
	rm -f ambimux

clean:
	rm -f *o
	rm -f ambimux
	cd $(HTSDIR) && make clean && cd ../

