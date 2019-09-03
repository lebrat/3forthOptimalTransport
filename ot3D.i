/* rt.i */

%module ot3D 


%{
#define SWIG_FILE_WITH_INIT
#include "ot3D.h"
#include "optim.h"
%}

%include "numpy.i"
%include "std_list.i"
%include "std_vector.i"

%apply (double* IN_ARRAY1, int DIM1) {(double* box, int n)}
%apply (double* IN_ARRAY1, int DIM1) {(double* positionX, int nX)}
%apply (double* IN_ARRAY1, int DIM1) {(double* positionY, int nY)}
%apply (double* IN_ARRAY1, int DIM1) {(double* positionZ, int nZ)}
%apply (double* IN_ARRAY1, int DIM1) {(double* weight, int nW)}

%apply (int DIM1, double* INPLACE_ARRAY1) {(int nMass,double* mass),(int nCost,double* cost),(int nxtDeriv,double* xtDeriv),(int nXt, double* xt),(int nBar,double* bar),(int nRho,double* rho),(int nStep,double* step),(int nScoring,double* scoring),(int lenAlpha,double* alpha),(int lenPsi,double* psi)};


%template(IntList) std::list<int>;
%template(DoubleList) std::list<double>;
%template(IntVector) std::vector<int>;
%template(adj) std::vector<std::list<int>>;
%template(EdgeList) std::list<std::vector<int>>;
%template(adjEdge) std::vector<std::list<std::vector<int>>>;
%template(ListDouble) std::vector<double>;
%template(dualPointList) std::list<std::vector<double>>;
%template(primalPtlist) std::list<PrimalPoint>;
%template(Linelist) std::list<Line> ;



%init %{
    import_array();
%}


%include "ot3D.h"
%include "optim.h"