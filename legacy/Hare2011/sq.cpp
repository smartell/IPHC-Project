// sq.cpp  -- return the square of objects.
//
# include <fvar.hpp>
//
double sq(_CONST double& x)
   {
   double xsq;
   xsq = x * x;
   return(xsq);
   }

dvariable sq(_CONST dvariable& x)
   {
  RETURN_ARRAYS_INCREMENT();
   dvariable xsq;
   xsq = x * x;
   RETURN_ARRAYS_DECREMENT();
   return(xsq);
   }

dvector sq(_CONST dvector& x)
   {
   int ibot = x.indexmin();
   int itop = x.indexmax();
   dvector xsq(ibot, itop);
   for (int i=ibot; i<=itop; i++)
      xsq(i) = x(i) * x(i);
   return(xsq);
   }

dvar_vector sq(_CONST dvar_vector& x)
   {
   RETURN_ARRAYS_INCREMENT();
   int ibot = x.indexmin();
   int itop = x.indexmax();
   dvar_vector xsq(ibot, itop);
   for (int i=ibot; i<=itop; i++)
      xsq(i) = x(i) * x(i);
   RETURN_ARRAYS_DECREMENT();
   return(xsq);
   }

dmatrix sq(_CONST dmatrix& x)
   {
   int ibot = x.rowmin();
   int itop = x.rowmax();
   int jbot = x.colmin();
   int jtop = x.colmax();
   dmatrix xsq(ibot, itop, jbot, jtop);
   for (int i=ibot; i<=itop; i++)
      for (int j=jbot; j<=jtop; j++)
         xsq(i,j) = x(i,j) * x(i,j);
   return(xsq);
   }

dvar_matrix sq(_CONST dvar_matrix& x)
   {
RETURN_ARRAYS_INCREMENT();
   int ibot = x.rowmin();
   int itop = x.rowmax();
   int jbot = x.colmin();
   int jtop = x.colmax();
   dvar_matrix xsq(ibot, itop, jbot, jtop);
   for (int i=ibot; i<=itop; i++)
      for (int j=jbot; j<=jtop; j++)
         xsq(i,j) = x(i,j) * x(i,j);
   RETURN_ARRAYS_DECREMENT();
   return(xsq);
   }
