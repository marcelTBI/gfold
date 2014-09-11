/* function from fold.c */
extern double  fold(const char *sequence, char *structure);
/* calculate mfe-structure of sequence */
extern double  energy_of_struct(const char *string, const char *structure);
/* calculate energy of string on structure */
extern void   free_arrays(void);           /* free arrays for mfe folding */
extern void   initialize_fold(int length); /* allocate arrays for folding */
extern void   update_fold_params(void);    /* recalculate parameters */
extern char  *backtrack_fold_from_pair(char *sequence, int i, int j);
extern double energy_of_circ_struct(const char *string, const char *structure);
extern int loop_energy(short * ptable, short *s, short *s1, int i);
