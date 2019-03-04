#include "include/rcore.h"
#include "include/sets.h"
#include "include/tests.h"
#include "include/dataframe.h"
#include "include/nmath.h"
#include "include/dpq.h"


// Define the maximum number of parent nodes.
#define MP 12


struct mi_xy
{
    int index_xy;
    int index_xy_x;
    int index_xy_y;
    int count;
}mi_xy[2000];

struct mi_x
{
    int index;
    int count;
}mi_x[2000];

struct mi_y
{
    int index;
    int count;
}mi_y[2000];

struct cmi_xyz
{
    double index_xyz;
    double index_xyz_xz;
    double index_xyz_yz;
    int index_xyz_z;
    int count;
}cmi_xyz[2000];

struct cmi_xz
{
    double index;
    int count;
}cmi_xz[2000];

struct cmi_yz
{
    double index;
    int count;
}cmi_yz[2000];

struct cmi_z
{
    int index;
    int count;
}cmi_z[2];


//  combination --> index --> count
//  SNPset: a subset of parent nodes' index

void GetMiCount(int *SNPset[], struct mi_xy *mi_xy, struct mi_x *mi_x, struct mi_y *mi_y, int k, int samplesize)
{
    int i = 0, j = 0, temp_xy = 0, temp_x = 0, temp_y = 0, n = 0, sign_xy = 0, sign_x = 0, sign_y = 0, sum_index = 0, temp_index = 0;

    for (i = 0; i < samplesize; i++)
    {
        // count for xy        
        sum_index = 0;
        for (j = 0; j < k; j++)
        {
            sum_index += (SNPset[j][i] - 1) * (int)pow(3.0, (k - 1 - j));
        }
        for (n = 0, sign_xy = 0; n < samplesize; n++)
        {
            if (sum_index == mi_xy[n].index_xy)
            {
                sign_xy = 1;
                temp_index = n;
                break;
            }
        }
        if (sign_xy == 0)
        {
            mi_xy[temp_xy].index_xy = sum_index;
            mi_xy[temp_xy].count++;
            temp_xy++;
        }
        else mi_xy[temp_index].count++;
        // count for x
        sum_index = 0;
        for (j = 0; j < MP; j++)
        {
            sum_index += (SNPset[j][i] - 1) * (int)pow(3.0, (MP - 1 - j));
        }
        for (n = 0, sign_x = 0; n < samplesize; n++)
        {
            if (sum_index == mi_x[n].index)
            {
                sign_x = 1;
                temp_index = n;
                break;
            }
        }
        if (sign_x == 0)
        {
            mi_x[temp_x].index = sum_index;
            mi_x[temp_x].count++;
            temp_x++;
            mi_xy[temp_xy - 1].index_xy_x = sum_index;
        }
        else 
        {
            mi_x[temp_index].count++;
            if (sign_xy == 0)
            {
                mi_xy[temp_xy - 1].index_xy_x = sum_index;
            }
        }
        // count for y
        sum_index = 0;
        for (j = 0; j < k - MP; j++)
        {
            sum_index += (SNPset[j + MP][i] - 1) * (int)pow(3.0, ((k - MP) - 1 - j));
        }        
        for (n = 0, sign_y = 0; n < samplesize; n++)
        {
            if (sum_index == mi_y[n].index)
            {
                sign_y = 1;
                temp_index = n;
                break;
            }
        }
        if (sign_y == 0)
        {
            mi_y[temp_y].index = sum_index;
            mi_y[temp_y].count++;
            temp_y++;
            mi_xy[temp_xy - 1].index_xy_y = sum_index;
        }
        else 
        {
            mi_y[temp_index].count++;
            if (sign_xy == 0)
            {
                mi_xy[temp_xy - 1].index_xy_y = sum_index;
            }
        }
    }
}


void GetCmiCount(int *SNPset[], struct cmi_xyz *cmi_xyz, struct cmi_xz *cmi_xz, struct cmi_yz *cmi_yz, struct cmi_z *cmi_z, int k, int samplesize)
{
    int i = 0, j = 0, n = 0, temp_xyz = 0, temp_xz = 0, temp_yz = 0, sign_xyz = 0, sign_xz = 0, sign_yz = 0, temp_index = 0, index_class = 0;
    double sum_index = 0.0;

    for (i = 0; i < samplesize; i++)
    {
        index_class = SNPset[k - 1][i] - 1;
        // count for xyz        
        sum_index = 0;
        for (j = 0; j < k - 1; j++)
        {
            sum_index += (SNPset[j][i] - 1) * (int)pow(3.0, ((k - 1) - 1 - j));
        }
        if (index_class == 1)
            sum_index += 0.5;

        for (n = 0, sign_xyz = 0; n < samplesize; n++)
        {
            if (sum_index == cmi_xyz[n].index_xyz)
            {
                sign_xyz = 1;
                temp_index = n;
                break;
            }
        }
        if (sign_xyz == 0)
        {
            cmi_xyz[temp_xyz].index_xyz = sum_index;
            cmi_xyz[temp_xyz].count++;
            temp_xyz++;
        }
        else cmi_xyz[temp_index].count++;
        // count for xz
        sum_index = 0;
        for (j = 0; j < MP; j++)
        {
            sum_index += (SNPset[j][i] - 1) * (int)pow(3.0, (MP - 1 - j));
        }
        if (index_class == 1)
            sum_index += 0.5;
        
        for (n = 0, sign_xz = 0; n < samplesize; n++)
        {
            if (sum_index == cmi_xz[n].index)
            {
                sign_xz = 1;
                temp_index = n;
                break;
            }
        }
        if (sign_xz == 0)
        {
            cmi_xz[temp_xz].index = sum_index;
            cmi_xz[temp_xz].count++;
            temp_xz++;
            cmi_xyz[temp_xyz - 1].index_xyz_xz = sum_index;
        }
        else 
        {
            cmi_xz[temp_index].count++;
            if (sign_xyz == 0)
            {
                cmi_xyz[temp_xyz - 1].index_xyz_xz = sum_index;
            }
        }
        // count for yz
        sum_index = 0;
        for (j = 0; j < k - MP - 1; j++)
        {
            sum_index += (SNPset[j + MP][i] - 1) * (int)pow(3.0, ((k - MP - 1) - 1 - j));
        } 
        if (index_class == 1)
            sum_index += 0.5;

        for (n = 0, sign_yz = 0; n < samplesize; n++)
        {
            if (sum_index == cmi_yz[n].index)
            {
                sign_yz = 1;
                temp_index = n;
                break;
            }
        }
        if (sign_yz == 0)
        {
            cmi_yz[temp_yz].index = sum_index;
            cmi_yz[temp_yz].count++;
            temp_yz++;
            cmi_xyz[temp_xyz - 1].index_xyz_yz = sum_index;
        }
        else
        {
            cmi_yz[temp_index].count++;
            if (sign_xyz == 0)
            {
                cmi_xyz[temp_xyz - 1].index_xyz_yz = sum_index;
            }
        }     
        // count for z
        if (index_class == 1)
        {
            cmi_z[1].count++;
            if (sign_xyz == 0)
                cmi_xyz[temp_xyz - 1].index_xyz_z = 1;
        }
        else 
        {
            cmi_z[0].count++;
            if (sign_xyz == 0)
                cmi_xyz[temp_xyz - 1].index_xyz_z = 0;
        }
    }
}


double dlik(SEXP x, double *nparams) {

int i = 0;
int *n = NULL, *xx = INTEGER(x), llx = NLEVELS(x), num = length(x);
double res = 0;

  /* initialize the contingency table. */
  fill_1d_table(xx, &n, llx, num);

  /* compute the entropy from the marginal frequencies. */
  for (i = 0; i < llx; i++)
    if (n[i] != 0)
      res += (double)n[i] * log((double)n[i] / num);

  /* we may want to store the number of parameters. */
  if (nparams)
    *nparams = llx - 1;

  Free1D(n);

  /* The node which do not have parents is so poor. */
  return res;

}/*DLIK*/

double mdlik(SEXP x, SEXP y, double *nparams) {

int i = 0, j = 0, k = 0;
int **n = NULL, *nj = NULL, *nk = NULL;
int llx = NLEVELS(x), lly = NLEVELS(y), num = length(x);
int *xx = INTEGER(x), *yy = INTEGER(y);
double res = 0;

  /* initialize the contingency table and the marginal frequencies. */
  n = (int **) Calloc2D(llx, lly, sizeof(int));
  nj = Calloc1D(lly, sizeof(int));
  nk = Calloc1D(llx, sizeof(int));
  
  /* compute the joint frequency of x and y. */
  /* caculate configurations of x and y. */
  for (k = 0; k < num; k++) 
    n[xx[k] - 1][yy[k] - 1]++;

  /* compute the marginals. */
  for (i = 0; i < llx; i++)
    for (j = 0; j < lly; j++)
      nj[j] += n[i][j];

  /* compute the marginals. */
  for (i = 0; i < llx; i++)
    for (j = 0; j < lly; j++)
      nk[i] += n[i][j];

  /* compute the conditional entropy from the joint and marginal
       frequencies. */
  for (i = 0; i < llx; i++)
    for (j = 0; j < lly; j++)
      if (n[i][j] != 0) 
        res += n[i][j] * log((double)n[i][j] * num / (nj[j] * nk[i]));
      
  /* we may want to store the number of parameters. */
  if (nparams)
    *nparams = (llx - 1) * lly;

  Free1D(nj);
  Free1D(nk);
  Free2D(n, llx);

  return res;

}/*CDLIK*/

double loglik_mit(SEXP target, SEXP x, SEXP data, double *nparams, int debuglevel) {

int i = 0, j = 0, df = 0;
double loglik = 0, pvalue = 0;
char *t = (char *)CHAR(STRING_ELT(target, 0));
SEXP nodes, node_t, parents, data_t, parent_vars, config;

/* get the node cached information. */
nodes = getListElement(x, "nodes");
node_t = getListElement(nodes, t);
/* get the parents of the node("Class" is always the last node). */
parents = getListElement(node_t, "parents");
/* extract the node's column from the data frame. */
data_t = c_dataframe_column(data, target, TRUE, FALSE);

if (length(parents) == 0)
{

  // loglik = dlik(data_t, nparams);
  // Rprintf("no parent, loglik = 0\n");
  loglik = 0.0;

  }/*THEN*/
  else {

    /* generate the configurations of the parents. */
    PROTECT(parent_vars = c_dataframe_column(data, parents, FALSE, FALSE));
    PROTECT(config = c_configurations(parent_vars, TRUE, TRUE));
    /* compute the log-likelihood. */
    loglik = mdlik(data_t, config, nparams);
    UNPROTECT(2);

    // Rprintf("Has %d parents, loglik = %lf\n", length(parents), loglik);

    for (i = 0; i < length(parents); i++) {

      if (i != 0) {

        // Rprintf("parents levels are %d\n", NLEVELS(VECTOR_ELT(parent_vars, i)));

        df = (NLEVELS(data_t) - 1) * (NLEVELS(VECTOR_ELT(parent_vars, i)) - 1);

        for (j = 0; j < i; j++)
          df *= NLEVELS(VECTOR_ELT(parent_vars, j));
        
      }/*THEN*/
      else 
        df = (NLEVELS(data_t) - 1) * (NLEVELS(config) - 1);

      pvalue += qgamma(0.95, 0.5 * df, 2.0, 1, 0);

    }/*FOR*/

    loglik = loglik * 2 - pvalue;

  }/*ELSE*/

  if (debuglevel > 0)
    Rprintf("  > loglikelihood is %lf.\n", loglik);

  return loglik;

}/*LOGLIK_DNODE*/

double cdlik(SEXP x, SEXP y, double *nparams) {

int i = 0, j = 0, k = 0;
int **n = NULL, *nj = NULL;
int llx = NLEVELS(x), lly = NLEVELS(y), num = length(x);
int *xx = INTEGER(x), *yy = INTEGER(y);
double res = 0;

  /* initialize the contingency table and the marginal frequencies. */
  n = (int **) Calloc2D(llx, lly, sizeof(int));
  nj = Calloc1D(lly, sizeof(int));

  /* compute the joint frequency of x and y. */
  for (k = 0; k < num; k++)
    n[xx[k] - 1][yy[k] - 1]++;

  /* compute the marginals. */
  for (i = 0; i < llx; i++)
    for (j = 0; j < lly; j++)
      nj[j] += n[i][j];

  /* compute the conditional entropy from the joint and marginal
       frequencies. */
  for (i = 0; i < llx; i++)
    for (j = 0; j < lly; j++)
      if (n[i][j] != 0)
        res += (double)n[i][j] * log((double)n[i][j] / (double)nj[j]);

  /* we may want to store the number of parameters. */
  if (nparams)
    *nparams = (llx - 1) * lly;

  Free1D(nj);
  Free2D(n, llx);

  return res;

}/*CDLIK*/



double loglik_dnode(SEXP target, SEXP x, SEXP data, double *nparams, int debuglevel) {

int i = 0, j = 0, llx = 0, lly = 0, llz = 0, length_xy = 0, samplesize = 0;
int *SNPdata_mxy[30] = {NULL}, *SNPdata_cmxyz[50] = {NULL}, *SNPdata_cmz[1] = {NULL};
double loglik = 0, loglik_temp = 0, nparams_temp = 1, inter = 0;
double sum_mi = 0, sum_cmi = 0, n_xy = 0, n_x = 0, n_y = 0, n_xyz = 0, n_xz = 0, n_yz = 0, n_z = 0;
char *t = (char *)CHAR(STRING_ELT(target, 0));
SEXP nodes, node_t, parents, data_t, parent_vars, parent_vars_MP, config, config_MP, parent_vars_all, parents_MP, parents_temp;


  /* get the node cached information. */
  nodes = getListElement(x, "nodes");
  node_t = getListElement(nodes, t);
  /* get the parents of the node. */
  parents = getListElement(node_t, "parents");
  /* extract the node's column from the data frame. */
  data_t = c_dataframe_column(data, target, TRUE, FALSE);
  SNPdata_cmz[0] = INTEGER(data_t); 
  samplesize = length(data_t);

  /* get the number of values of target. */
  llz = NLEVELS(data_t);


  if (length(parents) == 0) {

    loglik = dlik(data_t, nparams);

  }/*THEN*/
  else if (length(parents) <= MP && length(parents) > 0) {

    /* generate the configurations of the parents. */
    PROTECT(parent_vars = c_dataframe_column(data, parents, FALSE, FALSE));
    PROTECT(config = c_configurations(parent_vars, TRUE, TRUE));
    /* compute the log-likelihood. */
    loglik = cdlik(data_t, config, nparams);

    UNPROTECT(2);

  }/*ELSE*/
  else {
      do {
          /* generate the configurations of the parents. */
          PROTECT(parent_vars_all = c_dataframe_column(data, parents, FALSE, FALSE));

          length_xy = length(parents);

          for (i = 0; i < length(parents); i++)
          {
              SNPdata_mxy[i] = INTEGER(VECTOR_ELT(parent_vars_all, i));
              SNPdata_cmxyz[i] = INTEGER(VECTOR_ELT(parent_vars_all, i));
          }
          SNPdata_cmxyz[length(parents)] = SNPdata_cmz[0];

          /* we need to choose first MP parent nodes, and put them into parents_MP. */
          parents_MP = PROTECT(allocVector(STRSXP, MP));
          for (i = 0; i < MP; i++)
              SET_STRING_ELT(parents_MP, i, STRING_ELT(parents, i));
                            
          /* just a temporary variable. */
          parents_temp = PROTECT(allocVector(STRSXP, length(parents) - MP));
          for (i = MP, j = 0; i < length(parents); i++, j++)
              SET_STRING_ELT(parents_temp, j, STRING_ELT(parents, i));

          /* we need to update the parent nodes of parents. */
          parents = PROTECT(allocVector(STRSXP, length(parents_temp)));
          for (i = 0; i < length(parents_temp); i++)
              SET_STRING_ELT(parents, i, STRING_ELT(parents_temp, i));

          /* generate the configurations of the parents and the parents_MP. */
          PROTECT(parent_vars_MP = c_dataframe_column(data, parents_MP, FALSE, FALSE));
          PROTECT(config_MP = c_configurations(parent_vars_MP, TRUE, TRUE));
          PROTECT(parent_vars = c_dataframe_column(data, parents, FALSE, FALSE));
          PROTECT(config = c_configurations(parent_vars, TRUE, TRUE));
          /* compute the log-likelihood. */
          loglik = cdlik(data_t, config_MP, nparams);

          /* trying to caculate inter formula. */
          llx = NLEVELS(config_MP);
          lly = NLEVELS(config);
          inter = log(samplesize) * (llz - 1) * (llx + lly - llx * lly - 1) / 2 - dlik(data_t, nparams);

          for (i = 0; i < 2000; i++)
          {
              mi_xy[i].index_xy = -1;
              mi_xy[i].index_xy_x = -1;
              mi_xy[i].index_xy_y = -1;
              mi_xy[i].count = 0;
          }
          for (i = 0; i < 2000; i++)
          {
              cmi_xyz[i].index_xyz = -1;
              cmi_xyz[i].index_xyz_xz = -1;
              cmi_xyz[i].index_xyz_yz = -1;
              cmi_xyz[i].count = 0;
          }
          for (i = 0; i < 2000; i++)
          {
              mi_x[i].index = -1;
              mi_x[i].count = 0;
          }
          for (i = 0; i < 2000; i++)
          {
              mi_y[i].index = -1;
              mi_y[i].count = 0;
          }
          for (i = 0; i < 2000; i++)
          {
              cmi_xz[i].index = -1;
              cmi_xz[i].count = 0;
          }
          for (i = 0; i < 2000; i++)
          {
              cmi_yz[i].index = -1;
              cmi_yz[i].count = 0;
          }
          cmi_z[0].index = 0;
          cmi_z[0].count = 0;
          cmi_z[1].index = 1;
          cmi_z[1].count = 0;

          GetMiCount(SNPdata_mxy, mi_xy, mi_x, mi_y, length_xy, samplesize);
          GetCmiCount(SNPdata_cmxyz, cmi_xyz, cmi_xz, cmi_yz, cmi_z, length_xy + 1, samplesize);

          // mi calculator
          for (i = 0, sum_mi = 0; i < sizeof(mi_xy)/sizeof(struct mi_xy); i++)
          {
              if (mi_xy[i].index_xy != -1)
              {
                  n_xy = mi_xy[i].count;
                  for (j = 0; j < sizeof(mi_x)/sizeof(struct mi_x); j++)
                  {
                      if (mi_xy[i].index_xy_x == mi_x[j].index)
                        break;
                  }
                  n_x = mi_x[j].count;
                  for (j = 0; j < sizeof(mi_y)/sizeof(struct mi_y); j++)
                  {
                      if (mi_xy[i].index_xy_y == mi_y[j].index)
                        break;
                  }
                  n_y = mi_y[j].count;
                sum_mi += (n_xy * log(n_xy * samplesize / (n_x * n_y))) / samplesize;
              }
          }

          // cmi calculator
          for (i = 0, sum_cmi = 0; i < sizeof(cmi_xyz) / sizeof(struct cmi_xyz); i++)
          {
              if (cmi_xyz[i].index_xyz != -1)
              {
                  n_xyz = cmi_xyz[i].count;
                  for (j = 0; j < sizeof(cmi_xz) / sizeof(struct cmi_xz); j++)
                  {
                      if (cmi_xyz[i].index_xyz_xz == cmi_xz[j].index)
                          break;
                  }
                  n_xz = cmi_xz[j].count;
                  for (j = 0; j < sizeof(cmi_yz) / sizeof(struct cmi_yz); j++)
                  {
                      if (cmi_xyz[i].index_xyz_yz == cmi_yz[j].index)
                          break;
                  }
                  n_yz = cmi_yz[j].count;
                  for (j = 0; j < sizeof(cmi_z)/sizeof(struct cmi_z); j++)
                  {
                      if (cmi_xyz[i].index_xyz_z == cmi_z[j].index)
                        break;
                  }
                  n_z = cmi_z[j].count;
                  sum_cmi += (n_xyz * log(n_xyz * n_z / (n_xz * n_yz))) / samplesize;
              }
          }

          /* add them all together. */
          loglik_temp += loglik + inter + samplesize * (sum_cmi - sum_mi);        
          /* we don't use the default nparams, we caculate it by ourselves. */
          nparams_temp += llx;

          UNPROTECT(7);

      } while (length(parents) > MP);

      /* generate the configurations of the parents. */
      PROTECT(parent_vars = c_dataframe_column(data, parents, FALSE, FALSE));
      PROTECT(config = c_configurations(parent_vars, TRUE, TRUE));
      /* compute the log-likelihood. */
      loglik = cdlik(data_t, config, nparams);
      loglik_temp += loglik;
      lly = NLEVELS(config);
      nparams_temp += lly;
      *nparams = nparams_temp * (llz - 1);
      loglik = loglik_temp;

      UNPROTECT(2);

  }

  if (debuglevel > 0)
    Rprintf("  > loglikelihood is %lf.\n", loglik);

  return loglik;

}/*LOGLIK_DNODE*/




/*----------- DEBUGGING -------------
 * make CFLAGS='-DDEBUG_p -g'
 * (cd `R-devel RHOME`/src/nmath; gcc -I. -I../../src/include -I../../../R/src/include  -DHAVE_CONFIG_H -fopenmp -DDEBUG_p -g -c ../../../R/src/nmath/pgamma.c -o pgamma.o)
 */

/* Scalefactor:= (2^32)^8 = 2^256 = 1.157921e+77 */
#define SQR(x) ((x)*(x))
static const double scalefactor = SQR(SQR(SQR(4294967296.0)));
#undef SQR

/* If |x| > |k| * M_cutoff,  then  log[ exp(-x) * k^x ]	 =~=  -x */
static const double M_cutoff = M_LN2 * DBL_MAX_EXP / DBL_EPSILON;/*=3.196577e18*/

/* Continued fraction for calculation of
 *    1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3*d) + ... = sum_{k=0}^Inf x^k/(i+k*d)
 *
 * auxilary in log1pmx() and lgamma1p()
 */
static double
logcf (double x, double i, double d,
       double eps /* ~ relative tolerance */)
{
    double c1 = 2 * d;
    double c2 = i + d;
    double c4 = c2 + d;
    double a1 = c2;
    double b1 = i * (c2 - i * x);
    double b2 = d * d * x;
    double a2 = c4 * c2 - b2;

#if 0
    assert (i > 0);
    assert (d >= 0);
#endif

    b2 = c4 * b1 - i * b2;

    while (fabs(a2 * b1 - a1 * b2) > fabs(eps * b1 * b2)) {
	double c3 = c2*c2*x;
	c2 += d;
	c4 += d;
	a1 = c4 * a2 - c3 * a1;
	b1 = c4 * b2 - c3 * b1;

	c3 = c1 * c1 * x;
	c1 += d;
	c4 += d;
	a2 = c4 * a1 - c3 * a2;
	b2 = c4 * b1 - c3 * b2;

	if (fabs (b2) > scalefactor) {
	    a1 /= scalefactor;
	    b1 /= scalefactor;
	    a2 /= scalefactor;
	    b2 /= scalefactor;
	} else if (fabs (b2) < 1 / scalefactor) {
	    a1 *= scalefactor;
	    b1 *= scalefactor;
	    a2 *= scalefactor;
	    b2 *= scalefactor;
	}
    }

    return a2 / b2;
}

/* Accurate calculation of log(1+x)-x, particularly for small x.  */
double log1pmx (double x)
{
    static const double minLog1Value = -0.79149064;

    if (x > 1 || x < minLog1Value)
	return log1p(x) - x;
    else { /* -.791 <=  x <= 1  -- expand in  [x/(2+x)]^2 =: y :
	    * log(1+x) - x =  x/(2+x) * [ 2 * y * S(y) - x],  with
	    * ---------------------------------------------
	    * S(y) = 1/3 + y/5 + y^2/7 + ... = \sum_{k=0}^\infty  y^k / (2k + 3)
	   */
	double r = x / (2 + x), y = r * r;
	if (fabs(x) < 1e-2) {
	    static const double two = 2;
	    return r * ((((two / 9 * y + two / 7) * y + two / 5) * y +
			    two / 3) * y - x);
	} else {
	    static const double tol_logcf = 1e-14;
	    return r * (2 * y * logcf (y, 3, 2, tol_logcf) - x);
	}
    }
}


/* Compute  log(gamma(a+1))  accurately also for small a (0 < a < 0.5). */
double lgamma1p (double a)
{
    const double eulers_const =	 0.5772156649015328606065120900824024;

    /* coeffs[i] holds (zeta(i+2)-1)/(i+2) , i = 0:(N-1), N = 40 : */
    const int N = 40;
    static const double coeffs[40] = {
	0.3224670334241132182362075833230126e-0,/* = (zeta(2)-1)/2 */
	0.6735230105319809513324605383715000e-1,/* = (zeta(3)-1)/3 */
	0.2058080842778454787900092413529198e-1,
	0.7385551028673985266273097291406834e-2,
	0.2890510330741523285752988298486755e-2,
	0.1192753911703260977113935692828109e-2,
	0.5096695247430424223356548135815582e-3,
	0.2231547584535793797614188036013401e-3,
	0.9945751278180853371459589003190170e-4,
	0.4492623673813314170020750240635786e-4,
	0.2050721277567069155316650397830591e-4,
	0.9439488275268395903987425104415055e-5,
	0.4374866789907487804181793223952411e-5,
	0.2039215753801366236781900709670839e-5,
	0.9551412130407419832857179772951265e-6,
	0.4492469198764566043294290331193655e-6,
	0.2120718480555466586923135901077628e-6,
	0.1004322482396809960872083050053344e-6,
	0.4769810169363980565760193417246730e-7,
	0.2271109460894316491031998116062124e-7,
	0.1083865921489695409107491757968159e-7,
	0.5183475041970046655121248647057669e-8,
	0.2483674543802478317185008663991718e-8,
	0.1192140140586091207442548202774640e-8,
	0.5731367241678862013330194857961011e-9,
	0.2759522885124233145178149692816341e-9,
	0.1330476437424448948149715720858008e-9,
	0.6422964563838100022082448087644648e-10,
	0.3104424774732227276239215783404066e-10,
	0.1502138408075414217093301048780668e-10,
	0.7275974480239079662504549924814047e-11,
	0.3527742476575915083615072228655483e-11,
	0.1711991790559617908601084114443031e-11,
	0.8315385841420284819798357793954418e-12,
	0.4042200525289440065536008957032895e-12,
	0.1966475631096616490411045679010286e-12,
	0.9573630387838555763782200936508615e-13,
	0.4664076026428374224576492565974577e-13,
	0.2273736960065972320633279596737272e-13,
	0.1109139947083452201658320007192334e-13/* = (zeta(40+1)-1)/(40+1) */
    };

    const double c = 0.2273736845824652515226821577978691e-12;/* zeta(N+2)-1 */
    const double tol_logcf = 1e-14;
    double lgam;
    int i;

    if (fabs (a) >= 0.5)
	return lgammafn (a + 1);

    /* Abramowitz & Stegun 6.1.33 : for |x| < 2,
     * <==> log(gamma(1+x)) = -(log(1+x) - x) - gamma*x + x^2 * \sum_{n=0}^\infty c_n (-x)^n
     * where c_n := (Zeta(n+2) - 1)/(n+2)  = coeffs[n]
     *
     * Here, another convergence acceleration trick is used to compute
     * lgam(x) :=  sum_{n=0..Inf} c_n (-x)^n
     */
    lgam = c * logcf(-a / 2, N + 2, 1, tol_logcf);
    for (i = N - 1; i >= 0; i--)
	lgam = coeffs[i] - a * lgam;

    return (a * lgam - eulers_const) * a - log1pmx (a);
} /* lgamma1p */



/*
 * Compute the log of a sum from logs of terms, i.e.,
 *
 *     log (exp (logx) + exp (logy))
 *
 * without causing overflows and without throwing away large handfuls
 * of accuracy.
 */
double logspace_add (double logx, double logy)
{
    return fmax2 (logx, logy) + log1p (exp (-fabs (logx - logy)));
}


/*
 * Compute the log of a difference from logs of terms, i.e.,
 *
 *     log (exp (logx) - exp (logy))
 *
 * without causing overflows and without throwing away large handfuls
 * of accuracy.
 */
double logspace_sub (double logx, double logy)
{
    return logx + R_Log1_Exp(logy - logx);
}

/*
 * Compute the log of a sum from logs of terms, i.e.,
 *
 *     log (sum_i  exp (logx[i]) ) =
 *     log (e^M * sum_i  e^(logx[i] - M) ) =
 *     M + log( sum_i  e^(logx[i] - M)
 *
 * without causing overflows or throwing much accuracy.
 */
#ifdef HAVE_LONG_DOUBLE
# define EXP expl
# define LOG logl
#else
# define EXP exp
# define LOG log
#endif
double logspace_sum (const double* logx, int n)
{
    if(n == 0) return ML_NEGINF; // = log( sum(<empty>) )
    if(n == 1) return logx[0];
    if(n == 2) return logspace_add(logx[0], logx[1]);
    // else (n >= 3) :
    int i;
    // Mx := max_i log(x_i)
    double Mx = logx[0];
    for(i = 1; i < n; i++) if(Mx < logx[i]) Mx = logx[i];
    LDOUBLE s = (LDOUBLE) 0.;
    for(i = 0; i < n; i++) s += EXP(logx[i] - Mx);
    return Mx + (double) LOG(s);
}

/* dpois_wrap (x__1, lambda) := dpois(x__1 - 1, lambda);  where
 * dpois(k, L) := exp(-L) L^k / gamma(k+1)  {the usual Poisson probabilities}
 *
 * and  dpois*(.., give_log = TRUE) :=  log( dpois*(..) )
*/
static double
dpois_wrap (double x_plus_1, double lambda, int give_log)
{
#ifdef DEBUG_p
    REprintf (" dpois_wrap(x+1=%.14g, lambda=%.14g, log=%d)\n",
	      x_plus_1, lambda, give_log);
#endif
    if (!R_FINITE(lambda))
	return R_D__0;
    if (x_plus_1 > 1)
	return dpois_raw (x_plus_1 - 1, lambda, give_log);
    if (lambda > fabs(x_plus_1 - 1) * M_cutoff)
	return R_D_exp(-lambda - lgammafn(x_plus_1));
    else {
	double d = dpois_raw (x_plus_1, lambda, give_log);
#ifdef DEBUG_p
	REprintf ("  -> d=dpois_raw(..)=%.14g\n", d);
#endif
	return give_log
	    ? d + log (x_plus_1 / lambda)
	    : d * (x_plus_1 / lambda);
    }
}

/*
 * Abramowitz and Stegun 6.5.29 [right]
 */
static double
pgamma_smallx (double x, double alph, int lower_tail, int log_p)
{
    double sum = 0, c = alph, n = 0, term;

#ifdef DEBUG_p
    REprintf (" pg_smallx(x=%.12g, alph=%.12g): ", x, alph);
#endif

    /*
     * Relative to 6.5.29 all terms have been multiplied by alph
     * and the first, thus being 1, is omitted.
     */

    do {
	n++;
	c *= -x / n;
	term = c / (alph + n);
	sum += term;
    } while (fabs (term) > DBL_EPSILON * fabs (sum));

#ifdef DEBUG_p
    REprintf ("%5.0f terms --> conv.sum=%g;", n, sum);
#endif
    if (lower_tail) {
	double f1 = log_p ? log1p (sum) : 1 + sum;
	double f2;
	if (alph > 1) {
	    f2 = dpois_raw (alph, x, log_p);
	    f2 = log_p ? f2 + x : f2 * exp (x);
	} else if (log_p)
	    f2 = alph * log (x) - lgamma1p (alph);
	else
	    f2 = pow (x, alph) / exp (lgamma1p (alph));
#ifdef DEBUG_p
    REprintf (" (f1,f2)= (%g,%g)\n", f1,f2);
#endif
	return log_p ? f1 + f2 : f1 * f2;
    } else {
	double lf2 = alph * log (x) - lgamma1p (alph);
#ifdef DEBUG_p
	REprintf (" 1:%.14g  2:%.14g\n", alph * log (x), lgamma1p (alph));
	REprintf (" sum=%.14g  log(1+sum)=%.14g	 lf2=%.14g\n",
		  sum, log1p (sum), lf2);
#endif
	if (log_p)
	    return R_Log1_Exp (log1p (sum) + lf2);
	else {
	    double f1m1 = sum;
	    double f2m1 = expm1 (lf2);
	    return -(f1m1 + f2m1 + f1m1 * f2m1);
	}
    }
} /* pgamma_smallx() */

static double
pd_upper_series (double x, double y, int log_p)
{
    double term = x / y;
    double sum = term;

    do {
	y++;
	term *= x / y;
	sum += term;
    } while (term > sum * DBL_EPSILON);

    /* sum =  \sum_{n=1}^ oo  x^n     / (y*(y+1)*...*(y+n-1))
     *	   =  \sum_{n=0}^ oo  x^(n+1) / (y*(y+1)*...*(y+n))
     *	   =  x/y * (1 + \sum_{n=1}^oo	x^n / ((y+1)*...*(y+n)))
     *	   ~  x/y +  o(x/y)   {which happens when alph -> Inf}
     */
    return log_p ? log (sum) : sum;
}

/* Continued fraction for calculation of
 *    scaled upper-tail F_{gamma}
 *  ~=  (y / d) * [1 +  (1-y)/d +  O( ((1-y)/d)^2 ) ]
 */
static double
pd_lower_cf (double y, double d)
{
    double f= 0.0 /* -Wall */, of, f0;
    double i, c2, c3, c4,  a1, b1,  a2, b2;

#define	NEEDED_SCALE				\
	  (b2 > scalefactor) {			\
	    a1 /= scalefactor;			\
	    b1 /= scalefactor;			\
	    a2 /= scalefactor;			\
	    b2 /= scalefactor;			\
	}

#define max_it 200000

#ifdef DEBUG_p
    REprintf("pd_lower_cf(y=%.14g, d=%.14g)", y, d);
#endif
    if (y == 0) return 0;

    f0 = y/d;
    /* Needed, e.g. for  pgamma(10^c(100,295), shape= 1.1, log=TRUE): */
    if(fabs(y - 1) < fabs(d) * DBL_EPSILON) { /* includes y < d = Inf */
#ifdef DEBUG_p
	REprintf(" very small 'y' -> returning (y/d)\n");
#endif
	return (f0);
    }

    if(f0 > 1.) f0 = 1.;
    c2 = y;
    c4 = d; /* original (y,d), *not* potentially scaled ones!*/

    a1 = 0; b1 = 1;
    a2 = y; b2 = d;

    while NEEDED_SCALE

    i = 0; of = -1.; /* far away */
    while (i < max_it) {

	i++;	c2--;	c3 = i * c2;	c4 += 2;
	/* c2 = y - i,  c3 = i(y - i),  c4 = d + 2i,  for i odd */
	a1 = c4 * a2 + c3 * a1;
	b1 = c4 * b2 + c3 * b1;

	i++;	c2--;	c3 = i * c2;	c4 += 2;
	/* c2 = y - i,  c3 = i(y - i),  c4 = d + 2i,  for i even */
	a2 = c4 * a1 + c3 * a2;
	b2 = c4 * b1 + c3 * b2;

	if NEEDED_SCALE

	if (b2 != 0) {
	    f = a2 / b2;
	    /* convergence check: relative; "absolute" for very small f : */
	    if (fabs (f - of) <= DBL_EPSILON * fmax2(f0, fabs(f))) {
#ifdef DEBUG_p
		REprintf(" %g iter.\n", i);
#endif
		return f;
	    }
	    of = f;
	}
    }

    MATHLIB_WARNING(" ** NON-convergence in pgamma()'s pd_lower_cf() f= %g.\n",
		    f);
    return f;/* should not happen ... */
} /* pd_lower_cf() */
#undef NEEDED_SCALE


static double
pd_lower_series (double lambda, double y)
{
    double term = 1, sum = 0;

#ifdef DEBUG_p
    REprintf("pd_lower_series(lam=%.14g, y=%.14g) ...", lambda, y);
#endif
    while (y >= 1 && term > sum * DBL_EPSILON) {
	term *= y / lambda;
	sum += term;
	y--;
    }
    /* sum =  \sum_{n=0}^ oo  y*(y-1)*...*(y - n) / lambda^(n+1)
     *	   =  y/lambda * (1 + \sum_{n=1}^Inf  (y-1)*...*(y-n) / lambda^n)
     *	   ~  y/lambda + o(y/lambda)
     */
#ifdef DEBUG_p
    REprintf(" done: term=%g, sum=%g, y= %g\n", term, sum, y);
#endif

    if (y != floor (y)) {
	/*
	 * The series does not converge as the terms start getting
	 * bigger (besides flipping sign) for y < -lambda.
	 */
	double f;
#ifdef DEBUG_p
	REprintf(" y not int: add another term ");
#endif
	/* FIXME: in quite few cases, adding  term*f  has no effect (f too small)
	 *	  and is unnecessary e.g. for pgamma(4e12, 121.1) */
	f = pd_lower_cf (y, lambda + 1 - y);
#ifdef DEBUG_p
	REprintf("  (= %.14g) * term = %.14g to sum %g\n", f, term * f, sum);
#endif
	sum += term * f;
    }

    return sum;
} /* pd_lower_series() */

/*
 * Compute the following ratio with higher accuracy that would be had
 * from doing it directly.
 *
 *		 dnorm (x, 0, 1, FALSE)
 *	   ----------------------------------
 *	   pnorm (x, 0, 1, lower_tail, FALSE)
 *
 * Abramowitz & Stegun 26.2.12
 */
static double
dpnorm (double x, int lower_tail, double lp)
{
    /*
     * So as not to repeat a pnorm call, we expect
     *
     *	 lp == pnorm (x, 0, 1, lower_tail, TRUE)
     *
     * but use it only in the non-critical case where either x is small
     * or p==exp(lp) is close to 1.
     */

    if (x < 0) {
	x = -x;
	lower_tail = !lower_tail;
    }

    if (x > 10 && !lower_tail) {
	double term = 1 / x;
	double sum = term;
	double x2 = x * x;
	double i = 1;

	do {
	    term *= -i / x2;
	    sum += term;
	    i += 2;
	} while (fabs (term) > DBL_EPSILON * sum);

	return 1 / sum;
    } else {
	double d = dnorm (x, 0., 1., FALSE);
	return d / exp (lp);
    }
}

/*
 * Asymptotic expansion to calculate the probability that Poisson variate
 * has value <= x.
 * Various assertions about this are made (without proof) at
 * http://members.aol.com/iandjmsmith/PoissonApprox.htm
 */
static double
ppois_asymp (double x, double lambda, int lower_tail, int log_p)
{
    static const double coefs_a[8] = {
	-1e99, /* placeholder used for 1-indexing */
	2/3.,
	-4/135.,
	8/2835.,
	16/8505.,
	-8992/12629925.,
	-334144/492567075.,
	698752/1477701225.
    };

    static const double coefs_b[8] = {
	-1e99, /* placeholder */
	1/12.,
	1/288.,
	-139/51840.,
	-571/2488320.,
	163879/209018880.,
	5246819/75246796800.,
	-534703531/902961561600.
    };

    double elfb, elfb_term;
    double res12, res1_term, res1_ig, res2_term, res2_ig;
    double dfm, pt_, s2pt, f, np;
    int i;

    dfm = lambda - x;
    /* If lambda is large, the distribution is highly concentrated
       about lambda.  So representation error in x or lambda can lead
       to arbitrarily large values of pt_ and hence divergence of the
       coefficients of this approximation.
    */
    pt_ = - log1pmx (dfm / x);
    s2pt = sqrt (2 * x * pt_);
    if (dfm < 0) s2pt = -s2pt;

    res12 = 0;
    res1_ig = res1_term = sqrt (x);
    res2_ig = res2_term = s2pt;
    for (i = 1; i < 8; i++) {
	res12 += res1_ig * coefs_a[i];
	res12 += res2_ig * coefs_b[i];
	res1_term *= pt_ / i ;
	res2_term *= 2 * pt_ / (2 * i + 1);
	res1_ig = res1_ig / x + res1_term;
	res2_ig = res2_ig / x + res2_term;
    }

    elfb = x;
    elfb_term = 1;
    for (i = 1; i < 8; i++) {
	elfb += elfb_term * coefs_b[i];
	elfb_term /= x;
    }
    if (!lower_tail) elfb = -elfb;
#ifdef DEBUG_p
    REprintf ("res12 = %.14g   elfb=%.14g\n", elfb, res12);
#endif

    f = res12 / elfb;

    np = pnorm (s2pt, 0.0, 1.0, !lower_tail, log_p);

    if (log_p) {
	double n_d_over_p = dpnorm (s2pt, !lower_tail, np);
#ifdef DEBUG_p
	REprintf ("pp*_asymp(): f=%.14g	 np=e^%.14g  nd/np=%.14g  f*nd/np=%.14g\n",
		  f, np, n_d_over_p, f * n_d_over_p);
#endif
	return np + log1p (f * n_d_over_p);
    } else {
	double nd = dnorm (s2pt, 0., 1., log_p);

#ifdef DEBUG_p
	REprintf ("pp*_asymp(): f=%.14g	 np=%.14g  nd=%.14g  f*nd=%.14g\n",
		  f, np, nd, f * nd);
#endif
	return np + f * nd;
    }
} /* ppois_asymp() */


double pgamma_raw (double x, double alph, int lower_tail, int log_p)
{
/* Here, assume that  (x,alph) are not NA  &  alph > 0 . */

    double res;

#ifdef DEBUG_p
    REprintf("pgamma_raw(x=%.14g, alph=%.14g, low=%d, log=%d)\n",
	     x, alph, lower_tail, log_p);
#endif
    R_P_bounds_01(x, 0., ML_POSINF);

    if (x < 1) {
	res = pgamma_smallx (x, alph, lower_tail, log_p);
    } else if (x <= alph - 1 && x < 0.8 * (alph + 50)) {
	/* incl. large alph compared to x */
	double sum = pd_upper_series (x, alph, log_p);/* = x/alph + o(x/alph) */
	double d = dpois_wrap (alph, x, log_p);
#ifdef DEBUG_p
	REprintf(" alph 'large': sum=pd_upper*()= %.12g, d=dpois_w(*)= %.12g\n",
		 sum, d);
#endif
	if (!lower_tail)
	    res = log_p
		? R_Log1_Exp (d + sum)
		: 1 - d * sum;
	else
	    res = log_p ? sum + d : sum * d;
    } else if (alph - 1 < x && alph < 0.8 * (x + 50)) {
	/* incl. large x compared to alph */
	double sum;
	double d = dpois_wrap (alph, x, log_p);
#ifdef DEBUG_p
	REprintf(" x 'large': d=dpois_w(*)= %.14g ", d);
#endif
	if (alph < 1) {
	    if (x * DBL_EPSILON > 1 - alph)
		sum = R_D__1;
	    else {
		double f = pd_lower_cf (alph, x - (alph - 1)) * x / alph;
		/* = [alph/(x - alph+1) + o(alph/(x-alph+1))] * x/alph = 1 + o(1) */
		sum = log_p ? log (f) : f;
	    }
	} else {
	    sum = pd_lower_series (x, alph - 1);/* = (alph-1)/x + o((alph-1)/x) */
	    sum = log_p ? log1p (sum) : 1 + sum;
	}
#ifdef DEBUG_p
	REprintf(", sum= %.14g\n", sum);
#endif
	if (!lower_tail)
	    res = log_p ? sum + d : sum * d;
	else
	    res = log_p
		? R_Log1_Exp (d + sum)
		: 1 - d * sum;
    } else { /* x >= 1 and x fairly near alph. */
#ifdef DEBUG_p
	REprintf(" using ppois_asymp()\n");
#endif
	res = ppois_asymp (alph - 1, x, !lower_tail, log_p);
    }

    /*
     * We lose a fair amount of accuracy to underflow in the cases
     * where the final result is very close to DBL_MIN.	 In those
     * cases, simply redo via log space.
     */
    if (!log_p && res < DBL_MIN / DBL_EPSILON) {
	/* with(.Machine, double.xmin / double.eps) #|-> 1.002084e-292 */
#ifdef DEBUG_p
	REprintf(" very small res=%.14g; -> recompute via log\n", res);
#endif
	return exp (pgamma_raw (x, alph, lower_tail, 1));
    } else
	return res;
}


#ifdef DEBUG_qgamma
# define DEBUG_q
#endif

attribute_hidden
double qchisq_appr(double p, double nu, double g /* = log Gamma(nu/2) */,
		   int lower_tail, int log_p, double tol /* EPS1 */)
{
#define C7	4.67
#define C8	6.66
#define C9	6.73
#define C10	13.32

    double alpha, a, c, ch, p1;
    double p2, q, t, x;

    /* test arguments and initialise */

#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(nu))
	return p + nu;
#endif
    R_Q_P01_check(p);
    if (nu <= 0) ML_ERR_return_NAN;

    alpha = 0.5 * nu;/* = [pq]gamma() shape */
    c = alpha-1;

    if(nu < (-1.24)*(p1 = R_DT_log(p))) {	/* for small chi-squared */
	/* log(alpha) + g = log(alpha) + log(gamma(alpha)) =
	 *        = log(alpha*gamma(alpha)) = lgamma(alpha+1) suffers from
	 *  catastrophic cancellation when alpha << 1
	 */
	double lgam1pa = (alpha < 0.5) ? lgamma1p(alpha) : (log(alpha) + g);
	ch = exp((lgam1pa + p1)/alpha + M_LN2);
#ifdef DEBUG_qgamma
	REprintf(" small chi-sq., ch0 = %g\n", ch);
#endif

    } else if(nu > 0.32) {	/*  using Wilson and Hilferty estimate */

	x = qnorm(p, 0, 1, lower_tail, log_p);
	p1 = 2./(9*nu);
	ch = nu*pow(x*sqrt(p1) + 1-p1, 3);

#ifdef DEBUG_qgamma
	REprintf(" nu > .32: Wilson-Hilferty; x = %7g\n", x);
#endif
	/* approximation for p tending to 1: */
	if( ch > 2.2*nu + 6 )
	    ch = -2*(R_DT_Clog(p) - c*log(0.5*ch) + g);

    } else { /* "small nu" : 1.24*(-log(p)) <= nu <= 0.32 */

	ch = 0.4;
	a = R_DT_Clog(p) + g + c*M_LN2;
#ifdef DEBUG_qgamma
	REprintf(" nu <= .32: a = %7g\n", a);
#endif
	do {
	    q = ch;
	    p1 = 1. / (1+ch*(C7+ch));
	    p2 = ch*(C9+ch*(C8+ch));
	    t = -0.5 +(C7+2*ch)*p1 - (C9+ch*(C10+3*ch))/p2;
	    ch -= (1- exp(a+0.5*ch)*p2*p1)/t;
	} while(fabs(q - ch) > tol * fabs(ch));
    }

    return ch;
}

double qgamma(double p, double alpha, double scale, int lower_tail, int log_p)
/*			shape = alpha */
{
#define EPS1 1e-2
#define EPS2 5e-7/* final precision of AS 91 */
#define EPS_N 1e-15/* precision of Newton step / iterations */
#define LN_EPS -36.043653389117156 /* = log(.Machine$double.eps) iff IEEE_754 */

#define MAXIT 1000/* was 20 */

#define pMIN 1e-100   /* was 0.000002 = 2e-6 */
#define pMAX (1-1e-14)/* was (1-1e-12) and 0.999998 = 1 - 2e-6 */

    const static double
	i420  = 1./ 420.,
	i2520 = 1./ 2520.,
	i5040 = 1./ 5040;

    double p_, a, b, c, g, ch, ch0, p1;
    double p2, q, s1, s2, s3, s4, s5, s6, t, x;
    int i, max_it_Newton = 1;

    /* test arguments and initialise */

#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(alpha) || ISNAN(scale))
	return p + alpha + scale;
#endif
    R_Q_P01_boundaries(p, 0., ML_POSINF);

    if (alpha < 0 || scale <= 0) ML_ERR_return_NAN;

    if (alpha == 0) /* all mass at 0 : */ return 0.;

    if (alpha < 1e-10) {
    /* Warning seems unnecessary now: */
#ifdef _DO_WARN_qgamma_
	MATHLIB_WARNING(_("value of shape (%g) is extremely small: results may be unreliable"),
			alpha);
#endif
	max_it_Newton = 7;/* may still be increased below */
    }

    p_ = R_DT_qIv(p);/* lower_tail prob (in any case) */

#ifdef DEBUG_qgamma
    REprintf("qgamma(p=%7g, alpha=%7g, scale=%7g, l.t.=%2d, log_p=%2d): ",
	     p,alpha,scale, lower_tail, log_p);
#endif
    g = lgammafn(alpha);/* log Gamma(v/2) */

    /*----- Phase I : Starting Approximation */
    ch = qchisq_appr(p, /* nu= 'df' =  */ 2*alpha, /* lgamma(nu/2)= */ g,
		     lower_tail, log_p, /* tol= */ EPS1);
    if(!R_FINITE(ch)) {
	/* forget about all iterations! */
	max_it_Newton = 0; goto END;
    }
    if(ch < EPS2) {/* Corrected according to AS 91; MM, May 25, 1999 */
	max_it_Newton = 20;
	goto END;/* and do Newton steps */
    }

    /* FIXME: This (cutoff to {0, +Inf}) is far from optimal
     * -----  when log_p or !lower_tail, but NOT doing it can be even worse */
    if(p_ > pMAX || p_ < pMIN) {
	/* did return ML_POSINF or 0.;	much better: */
	max_it_Newton = 20;
	goto END;/* and do Newton steps */
    }

#ifdef DEBUG_qgamma
    REprintf("\t==> ch = %10g:", ch);
#endif

/*----- Phase II: Iteration
 *	Call pgamma() [AS 239]	and calculate seven term taylor series
 */
    c = alpha-1;
    s6 = (120+c*(346+127*c)) * i5040; /* used below, is "const" */

    ch0 = ch;/* save initial approx. */
    for(i=1; i <= MAXIT; i++ ) {
	q = ch;
	p1 = 0.5*ch;
	p2 = p_ - pgamma_raw(p1, alpha, /*lower_tail*/TRUE, /*log_p*/FALSE);
#ifdef DEBUG_qgamma
	if(i == 1) REprintf(" Ph.II iter; ch=%g, p2=%g\n", ch, p2);
	if(i >= 2) REprintf("     it=%d,  ch=%g, p2=%g\n", i, ch, p2);
#endif
#ifdef IEEE_754
	if(!R_FINITE(p2) || ch <= 0)
#else
	if(errno != 0 || ch <= 0)
#endif
	    { ch = ch0; max_it_Newton = 27; goto END; }/*was  return ML_NAN;*/

	t = p2*exp(alpha*M_LN2+g+p1-c*log(ch));
	b = t/ch;
	a = 0.5*t - b*c;
	s1 = (210+ a*(140+a*(105+a*(84+a*(70+60*a))))) * i420;
	s2 = (420+ a*(735+a*(966+a*(1141+1278*a)))) * i2520;
	s3 = (210+ a*(462+a*(707+932*a))) * i2520;
	s4 = (252+ a*(672+1182*a) + c*(294+a*(889+1740*a))) * i5040;
	s5 = (84+2264*a + c*(1175+606*a)) * i2520;

	ch += t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
	if(fabs(q - ch) < EPS2*ch)
	    goto END;
	if(fabs(q - ch) > 0.1*ch) {/* diverging? -- also forces ch > 0 */
	    if(ch < q) ch = 0.9 * q; else ch = 1.1 * q;
	}
    }
/* no convergence in MAXIT iterations -- but we add Newton now... */
#ifdef DEBUG_q
    MATHLIB_WARNING3("qgamma(%g) not converged in %d iterations; rel.ch=%g\n",
		     p, MAXIT, ch/fabs(q - ch));
#endif
/* was
 *    ML_ERROR(ME_PRECISION, "qgamma");
 * does nothing in R !*/

END:
/* PR# 2214 :	 From: Morten Welinder <terra@diku.dk>, Fri, 25 Oct 2002 16:50
   --------	 To: R-bugs@biostat.ku.dk     Subject: qgamma precision

   * With a final Newton step, double accuracy, e.g. for (p= 7e-4; nu= 0.9)
   *
   * Improved (MM): - only if rel.Err > EPS_N (= 1e-15);
   *		    - also for lower_tail = FALSE	 or log_p = TRUE
   *		    - optionally *iterate* Newton
   */
    x = 0.5*scale*ch;
    if(max_it_Newton) {
	/* always use log scale */
	if (!log_p) {
	    p = log(p);
	    log_p = TRUE;
	}
	if(x == 0) {
	    const double _1_p = 1. + 1e-7;
	    const double _1_m = 1. - 1e-7;
	    x = DBL_MIN;
	    p_ = pgamma(x, alpha, scale, lower_tail, log_p);
	    if(( lower_tail && p_ > p * _1_p) ||
	       (!lower_tail && p_ < p * _1_m))
		return(0.);
	    /* else:  continue, using x = DBL_MIN instead of  0  */
	}
	else
	    p_ = pgamma(x, alpha, scale, lower_tail, log_p);
	if(p_ == ML_NEGINF) return 0; /* PR#14710 */
	for(i = 1; i <= max_it_Newton; i++) {
	    p1 = p_ - p;
#ifdef DEBUG_qgamma
	    if(i == 1) REprintf("\n it=%d: p=%g, x = %g, p.=%g; p1=d{p}=%g\n",
				i, p, x, p_, p1);
	    if(i >= 2) REprintf("          x{it= %d} = %g, p.=%g, p1=d{p}=%g\n",
				i,    x, p_, p1);
#endif
	    if(fabs(p1) < fabs(EPS_N * p))
		break;
	    /* else */
	    if((g = dgamma(x, alpha, scale, log_p)) == R_D__0) {
#ifdef DEBUG_q
		if(i == 1) REprintf("no final Newton step because dgamma(*)== 0!\n");
#endif
		break;
	    }
	    /* else :
	     * delta x = f(x)/f'(x);
	     * if(log_p) f(x) := log P(x) - p; f'(x) = d/dx log P(x) = P' / P
	     * ==> f(x)/f'(x) = f*P / P' = f*exp(p_) / P' (since p_ = log P(x))
	     */
	    t = log_p ? p1*exp(p_ - g) : p1/g ;/* = "delta x" */
	    t = lower_tail ? x - t : x + t;
	    p_ = pgamma (t, alpha, scale, lower_tail, log_p);
	    if (fabs(p_ - p) > fabs(p1) ||
		(i > 1 && fabs(p_ - p) == fabs(p1)) /* <- against flip-flop */) {
		/* no improvement */
#ifdef DEBUG_qgamma
		if(i == 1 && max_it_Newton > 1)
		    REprintf("no Newton step done since delta{p} >= last delta\n");
#endif
		break;
	    } /* else : */
#ifdef Harmful_notably_if_max_it_Newton_is_1
	    /* control step length: this could have started at
	       the initial approximation */
	    if(t > 1.1*x) t = 1.1*x;
	    else if(t < 0.9*x) t = 0.9*x;
#endif
	    x = t;
	}
    }

    return x;
}