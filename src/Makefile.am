# edit:

# path to your distribution of ViennaRNA package (version 2.1.7 is tested to work)
VIENNARNA = $(HOME)/software/ViennaRNA-2.1.7

# do not edit:
AUTOMAKE_OPTIONS=foreign
bin_PROGRAMS=gfold
gfold_SOURCES=aln_util.h  energy_const.h  energy_par.c \
 energy_par.h  fold.c  fold.h  fold_vars.c \
 fold_vars.h  head.h  intloops.h  main.c \
 pair_mat.h  params.c  params.h  \
 utils.c  utils.h  gfold.c  \
 PS_dot.h pfunc.c \
 gfold_cmdline.c
gfold_SOURCES = head.h main.c \
 gfold.c  \
 pfunc.c \
 gfold_cmdline.c

BUILT_SOURCES = gfold_cmdline.c

main.o: gfold_cmdline.o 

gfold_cmdline.c: gfold.ggo
	gengetopt -i gfold.ggo

LDADD = $(OPENMP_CFLAGS) -lm $(VIENNARNA)/lib/libRNA.a
AM_CPPFLAGS = -I$(VIENNARNA)/H -I$(VIENNARNA)/lib
AM_CFLAGS = $(OPENMP_CFLAGS)
