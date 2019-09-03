#ifndef OPTIM_H
#define OPTIM_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/SparseCore>
#include <list>

typedef Eigen::Triplet<double> T;

struct matCOO{
    std::list<int> row;
    std::list<int> col;
    std::list<double> val;
    uint nr,nc;
    matCOO(){};
    matCOO(uint size){nr = size;nc=size;};
    matCOO(uint sizer,uint sizec){nr = sizer;nc=sizec;};
    matCOO(Eigen::SparseMatrix<double> const & m);
    Eigen::SparseMatrix<double> toSparseEigen();
    Eigen::SparseMatrix<double> Extract(std::vector<bool> Mask,std::vector<int> const & NewIndice, int size);
    void push(uint Prow,uint Pcol,double Pval){row.push_back(Prow);col.push_back(Pcol);val.push_back(Pval);};
    int getDimRow(){return (int) nr;};
    int getDimCol(){return (int) nc;};
    bool isSquared(){return nr == nc;};
};

class BFGS{
    const uint memorySize;
    const uint problemSize;
    std::list<Eigen::VectorXd> memoryX;
    std::list<Eigen::VectorXd> memoryGrad;
    public:
        BFGS();
        BFGS(uint memSize,uint pbSize);
        void addDirection(double step, Eigen::VectorXd direction, Eigen::VectorXd gradNew, Eigen::VectorXd gradOld);
        Eigen::VectorXd findDirection(Eigen::VectorXd gradf);

            
};

class gradientDescent{
    double stepSize;
    const double gradNormThres;
    const uint iterMax;
    public:
        gradientDescent();
        gradientDescent(double step, double gradN,int nbIter);
        void setStep(double s){ this->stepSize = s;};
        double step(){ return this->stepSize;};
        double gradN(){ return this->gradNormThres;};
        uint iMax(){ return this->iterMax;};

};

class gradientStrongWolf{
    const double stepSize;
    const double gradNormThres;
    const double alpha1;
    const double alpha2;
    const uint maxStepReduction;
    const uint iterMax;

    public:
        gradientStrongWolf();
        gradientStrongWolf(double step, double gradN,double wOne, double wTwo,
                           int innerLoop,int nbIter);
        
        double step(){ return this->stepSize;};
        double gradN(){ return this->gradNormThres;};
        double w1(){ return this->alpha1;};
        double w2(){ return this->alpha2;};
        uint innerMax(){ return this->maxStepReduction;};
        uint iMax(){ return this->iterMax;};

};

class barzilaiBorwein{
    const double gradNormThres;
    const uint iterMax;
    public:
        barzilaiBorwein();
        barzilaiBorwein(double gradN,int nbIter);
        double gradN(){ return this->gradNormThres;};
        uint iMax(){ return this->iterMax;};
        double bbStep;
};

class barzilaiBorweinLS{
    const double gradNormThres;
    const double lineSearchFactor;
    const double gammaWolfe;
    const uint iterMax;
    const uint memorySize;
    protected:
        std::list<double> costFunctionList;
    public:
        barzilaiBorweinLS();
        barzilaiBorweinLS(double gradN,double lineSearch,double gammaW,int nbIter,int memSize);
        double gradN(){ return this->gradNormThres;};
        double lSearch(){ return this->lineSearchFactor;};
        double gamma(){ return this->gammaWolfe;};
        void addCost(double cost);
        double minCost();
        uint iMax(){ return this->iterMax;};
        double bbStep;
};


class heavyBall{
    const double Lips;
    const double gradNormThres;
    const uint iterMax;
    public:
        heavyBall();
        heavyBall(double LispchitzConstant,double gradN,int nbIter);
        double L(){ return this->Lips;};
        double gradN(){ return this->gradNormThres;};
        uint iMax(){ return this->iterMax;};
};

#endif