AUTOMAKE_OPTIONS=foreign
bin_PROGRAMS=gfold
gfold_SOURCES=aln_util.h  energy_const.h  energy_par.c \
 energy_par.h  fold.c  fold.h  fold_vars.c \
 fold_vars.h  head.h  intloops.h  main.c \
 pair_mat.h  params.c  params.h  \
 utils.c  utils.h  gfold.c  \
 PS_dot.h pfunc.c \
 gfold_cmdline.c

BUILT_SOURCES = gfold_cmdline.c

gfold_cmdline.c: gfold.ggo
	gengetopt -i gfold.ggo

LDADD = -lm
