
CC = $(shell which gcc)

# CFLAGS = -g -O1 -Wall -Wextra -Wfloat-equal -Wno-unused-function -fsanitize=address -Wshadow -fno-omit-frame-pointer
CFLAGS = -g -O2 -Wall -Wextra -Wfloat-equal -Wno-unused-function -Wpointer-arith -Wshadow

HTSDIR = htslib
HTSLIB = $(HTSDIR)/libhts.a

CPPFLAGS += -I. -I$(HTSDIR)

all : ambimux

-include $(HTSDIR)/htslib.mk
-include $(HTSDIR)/htslib_static.mk

OBJS = main.o atac_data.o str_util.o sam_read.o \
	   variants.o bins.o overlap.o gtf_anno.o counts.o region.o \
	   clopts.o rna_data.o bc_stats.o mod.o bam_dat.o r_count.o \
	   bam_proc.o math_util.o array_util.o

LDFLAGS += -L$(HTSDIR)
LIBS += -lm -l:libhts.a -lpthread $(HTSLIB_static_LIBS)

ambimux_make : ambimux

hts :
	echo "building htslib"
	rm -rf $(HTSDIR)
	git clone --branch 1.17 https://github.com/samtools/htslib.git $(HTSDIR)
	cd $(HTSDIR) && git submodule update --init --recursive
	cd $(HTSDIR) && autoreconf -i && ./configure
	cd $(HTSDIR) && make lib-static

check_lib :
	if [ ! -d "$(HTSDIR)" ]; then \
		echo "htslib subdirectory not found, run 'make hts' first"; \
		exit 1; \
	fi

check_static :
	if [ ! -s "$(HTSDIR)/htslib_static.mk" ]; then \
		echo "htslib_static.mk not found, run 'make hts' first"; \
		exit 1; \
	fi

ambimux : check_lib check_static $(HTSDIR)/libhts.a $(OBJS)
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

