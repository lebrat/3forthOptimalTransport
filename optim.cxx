#include "optim.h"
#include <limits>
// #include <cstddef>
#include <boost/range/combine.hpp>
#include <boost/range/adaptor/reversed.hpp>



//////////////////////////////////////////////////////////////////
//////////////////////// COO Matrix Struct ///////////////////////
//////////////////////////////////////////////////////////////////

matCOO::matCOO(Eigen::SparseMatrix<double> const & m){
    nr = m.rows();
    nc = m.cols();
    for (int k=0; k < m.outerSize(); ++k){
        for (Eigen::SparseMatrix<double>::InnerIterator it(m,k); it; ++it)
        {
            row.push_back(it.row());
            col.push_back(it.col());
            val.push_back(it.value());
        }
    }
};

Eigen::SparseMatrix<double> matCOO::toSparseEigen(){
  std::vector<T> tripletList;
  tripletList.reserve(row.size());
  int i,j;
  double v;
  for(auto tup : boost::combine(row,col,val)){
    boost::tie(i,j,v) = tup;
    tripletList.push_back(T(i,j,v));
  }
  Eigen::SparseMatrix<double> mat(nr,nc);
  mat.setFromTriplets(tripletList.begin(), tripletList.end());
  return mat;  
};

Eigen::SparseMatrix<double> matCOO::Extract(std::vector<bool> Mask ,std::vector<int> const & NewIndice, int size){
    std::vector<T> tripletList;
    tripletList.reserve(row.size());
    int i,j;
    double v;
    for(auto tup : boost::combine(row,col,val)){   
        boost::tie(i,j,v) = tup;
        if (Mask[i] && Mask[j]){
            tripletList.push_back(T(NewIndice[i],NewIndice[j],v));
        }
        
    }
    Eigen::SparseMatrix<double> mat(size,size);
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return mat;  
};















//////////////////////////////////////////////////////////////////
/////////////////////////// BFGS CLASS ///////////////////////////
//////////////////////////////////////////////////////////////////

BFGS::BFGS() : memorySize(15), problemSize(1000){
}
BFGS::BFGS(uint memSize,uint pbSize) : memorySize(memSize), problemSize(pbSize){
    this->memoryX    = std::list<Eigen::VectorXd>();
    this->memoryGrad = std::list<Eigen::VectorXd>();
}
void BFGS::addDirection(double step, Eigen::VectorXd direction, Eigen::VectorXd gradNew, Eigen::VectorXd gradOld){
    if(this->memoryX.size()==memorySize)    this->memoryX.pop_back();
    if(this->memoryGrad.size()==memorySize) this->memoryGrad.pop_back();
    this->memoryX.push_front(step*direction);
    this->memoryGrad.push_front(gradOld - gradNew);
}

Eigen::VectorXd BFGS::findDirection(Eigen::VectorXd gradf){
    if(this->memoryX.size() == 0){
        return gradf;
    }
    else{
        Eigen::VectorXd r = -1.*gradf;
        std::vector<double> rhoVect = std::vector<double>(this->memoryX.size());
        Eigen::VectorXd x,grad;
        int compt = 0;
        double dotValue;
        for(auto tup : boost::combine(this->memoryX,this->memoryGrad)){
            boost::tie(x, grad) = tup;
            dotValue = grad.dot(x);
            rhoVect[compt] = std::abs(dotValue) < 1e-42 ? 0 : 1./dotValue;  
            compt++;
        }
        double gamma = this->memoryX.front().dot(this->memoryGrad.front())/(this->memoryGrad.front().dot(this->memoryGrad.front()));
        std::vector<double> alpha = std::vector<double>(this->memoryGrad.size());
        compt = 0; 
        for(auto tup : boost::combine(this->memoryX,this->memoryGrad)){
            boost::tie(x, grad) = tup;
            alpha[compt] = rhoVect[compt]*x.dot(r);
            r -= alpha[compt]*grad;
            compt++;
        }
        r *= gamma;
        double beta;
        compt = this->memoryX.size()-1;
        for(auto tup : boost::adaptors::reverse(boost::combine(this->memoryX,this->memoryGrad))){
            boost::tie(x, grad) = tup;
            beta = rhoVect[compt]*grad.dot(r);
            r += (alpha[compt] - beta)*x;
            compt--;
        }
        return -r;
        
    }			
}

//////////////////////////////////////////////////////////////////
////////////////// gradientDescent CLASS /////////////////////////
//////////////////////////////////////////////////////////////////
gradientDescent::gradientDescent() : stepSize(1e-2), gradNormThres(1e-5), iterMax(2000) {
}
gradientDescent::gradientDescent(double step, double gradN, int nbIter) : stepSize(step),
                                gradNormThres(gradN), iterMax((uint) nbIter){
}

//////////////////////////////////////////////////////////////////
///////////////// gradientStrongWolf CLASS ///////////////////////
//////////////////////////////////////////////////////////////////
gradientStrongWolf::gradientStrongWolf() : stepSize(1e-2), gradNormThres(1e-5),
                                           alpha1(0.001), alpha2(0.9), maxStepReduction(12), iterMax(2000) {
}
gradientStrongWolf::gradientStrongWolf(double step, double gradN, double wOne, double wTwo,int innerLoop, int nbIter) : 
                                stepSize(step), gradNormThres(gradN), alpha1(wOne), alpha2(wTwo),
                                maxStepReduction((uint) innerLoop), iterMax((uint) nbIter){
}

//////////////////////////////////////////////////////////////////
////////////////// barzilaiBorwein CLASS /////////////////////////
//////////////////////////////////////////////////////////////////
barzilaiBorwein::barzilaiBorwein() : gradNormThres(1e-5), iterMax(2000){
}
barzilaiBorwein::barzilaiBorwein(double gradN,int nbIter) : gradNormThres(gradN), iterMax((uint) nbIter) {
}

//////////////////////////////////////////////////////////////////
/////////////////// barzilaiBorweinLS CLASS //////////////////////
//////////////////////////////////////////////////////////////////
barzilaiBorweinLS::barzilaiBorweinLS() : gradNormThres(1e-5), lineSearchFactor(.70),gammaWolfe(1e-4), iterMax(2000),memorySize(10){
}
barzilaiBorweinLS::barzilaiBorweinLS(double gradN,double lineSearchFactor,double gammaW,int nbIter,int memSize) :
                 gradNormThres(gradN), lineSearchFactor(.70),gammaWolfe(gammaW),iterMax((uint) nbIter),memorySize((uint) memSize) {
}
void barzilaiBorweinLS::addCost(double cost){
    if(costFunctionList.size() > memorySize){
        costFunctionList.pop_back();
    }
    costFunctionList.push_front(cost);
}
double barzilaiBorweinLS::minCost(){
    if(costFunctionList.size() < 2){
        return std::numeric_limits<double>::min();
    }
    else{
        double minValue = -2.0*costFunctionList.front()-1.;
        for( const auto & elm : costFunctionList){
            if( elm < minValue ){
                minValue = elm;
            }
        }
        return minValue;
    }
}

//////////////////////////////////////////////////////////////////
/////////////////////// heavyBall CLASS //////////////////////////
//////////////////////////////////////////////////////////////////
heavyBall::heavyBall() : Lips(1000.), gradNormThres(1e-5), iterMax(2000){
}
heavyBall::heavyBall(double LispchitzConstant,double gradN,int nbIter) : Lips(LispchitzConstant),gradNormThres(gradN), iterMax((uint) nbIter) {
}
