/*****************************************************************************\
  This header file contains all the definitions and prototypes necessary for
 the routines in quad.c to solve the master problem using regularized QP method  which happens when config.MASTER_TYPE is equal to 1 (QP). 
 \*****************************************************************************/

#ifndef QUAD_H_
#define QUAD_H_

#include "sdconstants.h"
#include "sdglobal.h"

void change_bounds(prob_type *p, cell_type *c, soln_type *s);
void change_prob_type(one_problem *m);
void change_rhs_b(prob_type *p, cell_type *c, soln_type *s);
void change_rhs_back(prob_type *p, cell_type *c, soln_type *s);
void change_rhs_cuts(prob_type *p, cell_type *c, soln_type *s, cut_type *cuts);
void change_rhs(prob_type *p, cell_type *c, soln_type *s);
void construct_batch_QP(sdglobal_type* sd_global, prob_type *p, cell_type *c,
		double sigma);
void construct_QP(prob_type *p, cell_type *c, double sigma);

#endif
