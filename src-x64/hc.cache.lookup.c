#include "include/rcore.h"
#include "include/globals.h"
#include "include/graph.h"
#include "include/scores.h"
#include "include/matrix.h"
#include "include/learning.h"

SEXP score_cache_fill(SEXP nodes, SEXP data, SEXP network, SEXP score,
    SEXP extra, SEXP reference, SEXP equivalence, SEXP decomposability,
    SEXP updated, SEXP amat, SEXP cache, SEXP blmat, SEXP debug) {

int *colsum = NULL, nnodes = length(nodes), lupd = length(updated);
int *a = NULL, *upd = NULL, *b = NULL, debuglevel = isTRUE(debug);
int i = 0, j = 0, k = 0;
double *cache_value = NULL;
SEXP arc, delta, op, temp;

  /* save a pointer to the adjacency matrix, the blacklist and the
   * updated nodes. */
  a = INTEGER(amat);
  b = INTEGER(blmat);
  upd = INTEGER(updated);

  /* if there are no nodes to update, return. */
  if (lupd == 0) return cache;

  /* set up row and column total to check for score equivalence;
   * zero means no parent nodes. */
  if (isTRUE(equivalence)) {

    colsum = Calloc1D(nnodes, sizeof(int));

    for (i = 0; i < nnodes; i++)
      for (j = 0; j < nnodes; j++)
        colsum[j] += a[CMC(i, j, nnodes)];

  }/*THEN*/

  /* allocate and initialize the cache. */
  cache_value = REAL(cache);

  /* allocate a two-slot character vector. */
  PROTECT(arc = allocVector(STRSXP, 2));

  /* allocate and initialize the fake score delta. */
  PROTECT(delta = ScalarReal(0));

  /* allocate and initialize the score.delta() operator. */
  PROTECT(op = mkString("set"));

  for (i = 0; i < nnodes; i++) {

    for (j = 0; j < nnodes; j++) {

       /* incident nodes must be different from each other. */
       if (i == j) continue;

       /* if only one or two nodes' caches need updating, skip the rest. */
       for (k = 0; k < lupd; k++)
         if (upd[k] == j)
           goto there;

       continue;

there:

       /* no need to compute the score delta for blacklisted arcs. */
       if (b[CMC(i, j, nnodes)] == 1)
         continue;

       /* use score equivalence if possible to check only one orientation. */
       if (isTRUE(equivalence)) {

         /* if the following conditions are met, look up the score delta of
          * the reverse of the current arc:
          *   1) that score delta has already been computed.
          *   2) both incident nodes have no parent, so the arc is really
          *      score equivalent (no v-structures).
          *   3) the reversed arc has not been blacklisted, as the score delta
          *      is not computed in this case. */
         if ((i > j) && (colsum[i] + colsum[j] == 0) && (b[CMC(j, i, nnodes)] == 0)) {

           cache_value[CMC(i, j, nnodes)] = cache_value[CMC(j, i, nnodes)];
           continue;

         }/*THEN*/

       }/*THEN*/

       /* save the nodes incident on the arc. */
       SET_STRING_ELT(arc, 0, STRING_ELT(nodes, i));
       SET_STRING_ELT(arc, 1, STRING_ELT(nodes, j));

       /* if the arc is not present in the graph it should be added;
        * otherwise it should be removed. */
       if (a[CMC(i, j, nnodes)] == 0)
         SET_STRING_ELT(op, 0, mkChar("set"));
       else
         SET_STRING_ELT(op, 0, mkChar("drop"));

       /* checkpoint allocated memory. */
       /* evaluate the call to score.delta() for the arc. */
       PROTECT(temp = score_delta(arc, network, data, score, delta, reference,
         op, extra, decomposability));

       cache_value[CMC(i, j, nnodes)] = NUM(VECTOR_ELT(temp, 1));
       UNPROTECT(1);

       if (debuglevel > 0)
         Rprintf("* caching score delta for arc %s -> %s (%lf).\n",
           CHAR(STRING_ELT(nodes, i)), CHAR(STRING_ELT(nodes, j)),
            cache_value[CMC(i, j, nnodes)]);

    }/*FOR*/

  }/*FOR*/

  UNPROTECT(3);

  if (isTRUE(equivalence))
    Free1D(colsum);

  return cache;

}/*HC_CACHE_FILL*/

/* a single step of the optimized hill climbing (one arc addition/removal/reversal). */
SEXP hc_opt_step(SEXP amat, SEXP nodes, SEXP added, SEXP cache_bic, SEXP cache_mit,
                 SEXP reference_bic, SEXP reference_mit, SEXP wlmat, SEXP blmat,
                 SEXP nparents, SEXP maxp, SEXP debug) {

  int nnodes = length(nodes), i = 0, j = 0, opt_count = 0, index = 0;
  int *am = NULL, *ad = NULL, *w = NULL, *b = NULL, debuglevel = isTRUE(debug);
  int counter = 0, update = 1, from = 0, to = 0, *path = NULL, *scratch = NULL;
  double *cache_value_bic = NULL, *cache_value_mit = NULL, temp_bic = 0, temp_mit = 0, max = 0, tol = MACHINE_TOL;
  double *mp = REAL(maxp), *np = REAL(nparents);
  SEXP bestop;
  struct operation
  {
    int opt;
    char fr_node;
    char to_node;
    double up_score;
    double fitscore;
  }tempstore[100];

  /* allocate and initialize the return value (use FALSE as a canary value). */
  PROTECT(bestop = allocVector(VECSXP, 3));
  setAttrib(bestop, R_NamesSymbol, mkStringVec(3, "op", "from", "to"));

  /* allocate and initialize a dummy FALSE object. */
  SET_VECTOR_ELT(bestop, 0, ScalarLogical(FALSE));

  /* allocate buffers for c_has_path(). */
  path = Calloc1D(nnodes, sizeof(int));
  scratch = Calloc1D(nnodes, sizeof(int));

  /* save pointers to the numeric/integer matrices. */
  cache_value_bic = REAL(cache_bic);
  cache_value_mit = REAL(cache_mit);
  ad = INTEGER(added);
  am = INTEGER(amat);
  w = INTEGER(wlmat);
  b = INTEGER(blmat);

  if (debuglevel > 0) {

     /* count how may arcs are to be tested. */
     for (i = 0; i < nnodes * nnodes; i++)
       counter += ad[i];

     Rprintf("----------------------------------------------------------------\n");
     Rprintf("* trying to add one of %d arcs.\n", counter);

  }/*THEN*/

  for (i = 0; i < nnodes; i++) {

    for (j = 0; j < nnodes; j++) {

      /* nothing to see, move along. */
      if (ad[CMC(i, j, nnodes)] == 0)
        continue;

      /* retrieve the score delta from the cache. */
      temp_bic = cache_value_bic[CMC(i, j, nnodes)];
      temp_mit = cache_value_mit[CMC(i, j, nnodes)];

      if (debuglevel > 0) {

        Rprintf("  > trying to add %s -> %s.\n", NODE(i), NODE(j));
        Rprintf("    > delta between BIC scores for nodes %s %s is %lf.\n",
          NODE(i), NODE(j), temp_bic);
        Rprintf("    > delta between MIT scores for nodes %s %s is %lf.\n",
          NODE(i), NODE(j), temp_mit);          

      }/*THEN*/

      /* this score delta is the best one at the moment, so add the arc if it
       * does not introduce cycles in the graph. */
      if (temp_bic > 0 && temp_mit > 0)
      {

        if (c_has_path(j, i, am, nnodes, nodes, FALSE, FALSE, path, scratch,
                       FALSE))
        {

          if (debuglevel > 0)
            Rprintf("    > not adding, introduces cycles in the graph.\n");

          continue;

        } /*THEN*/

        max = temp_bic + temp_mit;
        tempstore[opt_count].opt = 1;
        tempstore[opt_count].fr_node = i;
        tempstore[opt_count].to_node = j;
        tempstore[opt_count].up_score = max;
        opt_count++;

      }/*THEN*/

    }/*FOR*/

  }/*FOR*/

  if (debuglevel > 0) {

     /* count how may arcs are to be tested. */
     for (i = 0, counter = 0; i < nnodes * nnodes; i++)
       counter += am[i] * (1 - w[i]);

     Rprintf("----------------------------------------------------------------\n");
     Rprintf("* trying to remove one of %d arcs.\n", counter);

  }/*THEN*/

  for (i = 0; i < nnodes; i++) {

    for (j = 0; j < nnodes; j++) {

      /* nothing to see, move along. */
      if (am[CMC(i, j, nnodes)] == 0)
        continue;

      /* whitelisted arcs are not to be removed, ever. */
      if (w[CMC(i, j, nnodes)] == 1)
        continue;

      /* retrieve the score delta from the cache. */
      temp_bic = cache_value_bic[CMC(i, j, nnodes)];
      temp_mit = cache_value_mit[CMC(i, j, nnodes)];

      if (debuglevel > 0) {

        Rprintf("  > trying to remove %s -> %s.\n", NODE(i), NODE(j));
        Rprintf("    > delta between BIC scores for nodes %s %s is %lf.\n",
          NODE(i), NODE(j), temp_bic);
        Rprintf("    > delta between MIT scores for nodes %s %s is %lf.\n",
          NODE(i), NODE(j), temp_mit);

      }/*THEN*/

      if (temp_bic > 0 && temp_mit > 0)
      {

        max = temp_bic + temp_mit;
        tempstore[opt_count].opt = 2;
        tempstore[opt_count].fr_node = i;
        tempstore[opt_count].to_node = j;
        tempstore[opt_count].up_score = max;
        opt_count++;

      }/*THEN*/

    }/*FOR*/

  }/*FOR*/

  if (debuglevel > 0) {

     /* count how may arcs are to be tested. */
     for (i = 0, counter = 0; i < nnodes; i++)
       for (j = 0; j < nnodes; j++)
         counter += am[CMC(i, j, nnodes)] * (1 - b[CMC(j, i, nnodes)]);

     Rprintf("----------------------------------------------------------------\n");
     Rprintf("* trying to reverse one of %d arcs.\n", counter);

  }/*THEN*/

  for (i = 0; i < nnodes; i++) {

    for (j = 0; j < nnodes; j++) {

      /* nothing to see, move along. */
      if (am[CMC(i, j, nnodes)] == 0)
        continue;

      /* don't reverse an arc if the one in the opposite direction is
       * blacklisted, ever. */
      if (b[CMC(j, i, nnodes)] == 1)
        continue;

      /* do not reverse an arc if that means violating the limit on the
       * maximum number of parents. */
      if (np[i] >= *mp)
        continue;

      /* retrieve the score delta from the cache. */
      temp_bic = cache_value_bic[CMC(i, j, nnodes)] + cache_value_bic[CMC(j, i, nnodes)];
      temp_mit = cache_value_mit[CMC(i, j, nnodes)] + cache_value_mit[CMC(j, i, nnodes)];
      /* nuke small values and negative zeroes. */
      if (fabs(temp_bic) < tol) temp_bic = 0;
      if (fabs(temp_mit) < tol) temp_mit = 0;

      if (debuglevel > 0) {

        Rprintf("  > trying to reverse %s -> %s.\n", NODE(i), NODE(j));
        Rprintf("    > delta between BIC scores for nodes %s %s is %lf.\n",
          NODE(i), NODE(j), temp_bic);
        Rprintf("    > delta between MIT scores for nodes %s %s is %lf.\n",
          NODE(i), NODE(j), temp_mit);

      }/*THEN*/

      if (temp_bic > 0 && temp_mit > 0)
      {

        if (c_has_path(i, j, am, nnodes, nodes, FALSE, TRUE, path, scratch,
                       FALSE))
        {

          if (debuglevel > 0)
            Rprintf("    > not reversing, introduces cycles in the graph.\n");

          continue;

        } /*THEN*/

        max = temp_bic + temp_mit;
        tempstore[opt_count].opt = 3;
        tempstore[opt_count].fr_node = i;
        tempstore[opt_count].to_node = j;
        tempstore[opt_count].up_score = max;
        opt_count++;

      }/*THEN*/

    }/*FOR*/

  }/*FOR*/

  for (i = 0, max = 0; i < opt_count; i++)
    max += tempstore[i].up_score;

  for (i = 0; i < opt_count; i++)
    tempstore[i].fitscore = tempstore[i].up_score / max; 

  double fitselect = rand()%100/(double)101;
  double fitsum = 0.0;

  for (i = 0; i < opt_count; i++)
  {
    fitsum += tempstore[i].fitscore;
    if (fitsum > fitselect)
    {
      index = i;
      break;
    }
  }

  from = tempstore[index].fr_node;
  to = tempstore[index].to_node;

  if (tempstore[index].opt == 1)
  {
    if (debuglevel > 0)
      Rprintf("    @ adding %s -> %s.\n", NODE(from), NODE(to));
    bestop_update(bestop, "set", NODE(from), NODE(to));
    update = 1;
  }
  if (tempstore[index].opt == 2)
  {
    if (debuglevel > 0)
      Rprintf("    @ removing %s -> %s.\n", NODE(from), NODE(to));
    bestop_update(bestop, "drop", NODE(from), NODE(to));
    update = 1;
  }
  if (tempstore[index].opt == 3)
  {
    if (debuglevel > 0)
      Rprintf("    @ reversing %s -> %s.\n", NODE(from), NODE(to));
    bestop_update(bestop, "reverse", NODE(from), NODE(to));
    update = 2;
  }

  /* update the reference scores. */
  REAL(reference_bic)[to] += cache_value_bic[CMC(from, to, nnodes)];
  if (update == 2)
    REAL(reference_bic)[from] += cache_value_bic[CMC(to, from, nnodes)];

  REAL(reference_mit)[to] += cache_value_mit[CMC(from, to, nnodes)];
  if (update == 2)
    REAL(reference_mit)[from] += cache_value_mit[CMC(to, from, nnodes)];

  Free1D(path);
  Free1D(scratch);

  UNPROTECT(1);

  return bestop;

}/*HC_OPT_STEP*/

void bestop_update(SEXP bestop, char *op, const char *from, const char *to) {

  SET_VECTOR_ELT(bestop, 0, mkString(op));
  SET_VECTOR_ELT(bestop, 1, mkString(from));
  SET_VECTOR_ELT(bestop, 2, mkString(to));

}/*BESTOP_UPDATE*/

SEXP hc_to_be_added(SEXP arcs, SEXP blacklist, SEXP whitelist, SEXP nparents,
    SEXP maxp, SEXP nodes, SEXP convert) {

int i = 0, j = 0, narcs = 0, dims = length(nodes);
int *a = NULL, *coords = NULL;
double *mp = REAL(maxp), *np = NULL;
short int referenced = 0;
SEXP try, result = R_NilValue, result2;

  /* transform the arc set into an adjacency matrix, if it's not one already. */
  if (isInteger(arcs)) {

    if ((referenced = MAYBE_REFERENCED(arcs)))
      PROTECT(result = duplicate(arcs));

  }/*THEN*/
  else {

    PROTECT(result = arcs2amat(arcs, nodes));

  }/*ELSE*/

  /* dereference the adjacency matrix once and for all. */
  a = INTEGER(result);

  /* compute the number the parents of each node, unless provided. */
  if (nparents == R_NilValue) {

    np = Calloc1D(dims, sizeof(double));
    for (i = 0; i < dims; i++)
      for (j = 0; j < dims; j++)
        np[j] = a[CMC(i, j, dims)];

  }/*THEN*/
  else {

    np = REAL(nparents);

  }/*ELSE*/

  /* flip all the nondiagonal cells. */
  for (j = 0; j < dims; j++) {

    for (i = 0; i < dims; i++) {

      /* diagonal elements are always equal to zero, skip them. */
      if (i == j)
        continue;

      a[CMC(i, j, dims)] = 1 - a[CMC(i, j, dims)];

    }/*FOR*/

  }/*FOR*/

  /* if an arc is present in the graph in one direction, you cannot add it in
   * the other direction (it would be a reversal); flip both in the adjacency
   * matrix. */
  for (j = 0; j < dims; j++)
    for (i = j + 1; i < dims; i++)
      a[CMC(j, i, dims)] = a[CMC(i, j, dims)] = a[CMC(i, j, dims)] * a[CMC(j, i, dims)];

  /* if a node has already reached its maximum number parents, do not add
   * more arcs pointing to that node. */
  for (j = 0; j < dims; j++)
    if (np[j] >= *mp)
      memset(a + j * dims, '\0', dims * sizeof(int));

#define FLIP_FROM_LIST(list, value) \
  if (!isNull(list)) { \
    if (!isInteger(list)) { \
      PROTECT(try = match(nodes, list, 0)); \
      coords = INTEGER(try); \
      narcs = length(try)/2; \
      for (i = 0; i < narcs; i++)  \
        a[CMC(coords[i] - 1, coords[i + narcs] - 1, dims)] = value; \
      UNPROTECT(1); \
    }/*THEN*/ \
    else { \
      coords = INTEGER(list); \
      for (i = 0; i < dims * dims; i ++) \
        if (coords[i] == 1) \
          a[i] = value; \
    }/*ELSE*/ \
  }/*THEN*/

  /* now the blacklist gets involved. */
  FLIP_FROM_LIST(blacklist, 0);
  /* and, last but not least, the whitelist gets involved. */
  FLIP_FROM_LIST(whitelist, 1);

  if (nparents == R_NilValue)
    Free1D(np);

  /* return either the adjacency matrix or the arc set. */
  if (isTRUE(convert)) {

    PROTECT(result2 = amat2arcs(result, nodes));

    if (referenced || !isInteger(arcs))
      UNPROTECT(2);
    else
      UNPROTECT(1);
    return result2;

  }/*THEN*/
  else {

    if (referenced || !isInteger(arcs))
      UNPROTECT(1);
    return result;

  }/*ELSE*/

}/*HC_TO_BE_ADDED*/

