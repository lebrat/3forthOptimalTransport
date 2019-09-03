#ifndef OT3D_H
#define OT3D_H


#include "optim.h"
#include "LaguerreCell3D.h"
#include <string>


struct derivInfo{
    Eigen::VectorXd gradient;
    double costFunction;
    int nbHidden;
    std::string ezHidden;
    derivInfo(){};
};

struct optimResult{
    Eigen::VectorXd phiOpt;
    double costFunction;
    int nbIter;
    double wallTime;
    optimResult(){};
};


struct optimResultResearch{
    std::list<double> stp;
    std::list<double> nrm;
    std::list<double> cst;
    std::list<double> tme;
    std::list<int> hddn;
    int nbtt;
    double wllTme;
    optimResultResearch(){};
};

class OTproblem{
    LaguerreCell3D Lag;
    Polyline pol;  
    int nPoints;
    derivInfo oracle(Eigen::VectorXd X);
    Eigen::VectorXd massDiracs;
    const std::string red,green,yellow,mage,cyan,blue,reset;
    optimResult optimizedPsi;
    optimResultResearch dataConvergence;
    public:
        bool parallelComputation;
        int nbCore;
        void printCGALVERSION();
        OTproblem():red("\033[0;31m"),green("\033[1;32m"),yellow("\033[1;33m"),mage("\033[1;35m"),
                    cyan("\033[1;36m"),blue("\033[0;34m"),reset("\033[0m"){};
        // position des points.
        void setLaguerrePosition(double* positionX,int nX,double* positionY,int nY,double* positionZ,int nZ);
        // polyLine avec sa masse.
        void setPolyline(double* positionX,int nX,double* positionY,int nY,double* positionZ,int nZ,double* weight,int nW);
        // Permet de construire la triangulation de laguerre assiociee
        // le weight ici est le parametre qui permet de faire bouger les 
        //  frontieres des cellules.
        void setPolylineNew(double* positionX,int nX,double* positionY,int nY,double* positionZ,int nZ,double* weight,int nW);

        void setLaguerrePsi(double* weight,int nW);
        // realise le calcul de la masse des couts locaux et des Barycentres.
        void computeAll(int nMass,double* mass,int nCost,double* cost,int nBar,double* bar);
        // calcul des derivees selon les points de la polyline.
        void computeDeriv(int nMass,double* mass,int nCost,double* cost);
        

        matCOO computeHessian();
        void computeDerivTotal(int nMass,double* mass,int nCost,double* cost,int nRho,double* rho);
        // calcul des infos de la polyline.
        void computePolyInfo(int nMass,double* mass,int nCost,double* cost);

        // some setters
        void setMassDiracs(int nMass,double* mass);

        // optimization methods
        void perfomOptimPsiBFGS(double gradTol,int nMaxIter,double w1, double w2, int wMaxIter,int memSize,bool parallel,int nbCore);
        void perfomOptimPsiConstantGD(double stepSize,double gradTol,int nMaxIter,bool parallel,int nbCore);
        void perfomOptimPsiGreedyGD(int nStep,double* step,int nScoring,double* scoring,double gradTol,int nMaxIter,bool parallel,int nbCore);
        void perfomOptimPsiStrongWolfGD(double stepSize,double gradTol,int nMaxIter,double w1, double w2, int wMaxIter,bool parallel,int nbCore);
        void perfomOptimPsiBarzilaiB(double gradTol,int nMaxIter,bool parallel,int nbCore);
        void perfomOptimPsiBarzilaiBWolf(double gradTol,double lineSearchRed,double gammaWolf,int nMaxIter,int memSize,bool parallel,int nbCore);
        void perfomOptimPsiBarzilaiBWolf2(double gradTol,double lineSearchRed,double gammaWolf,int nMaxIter,int memSize,int maxHiddenIncrease,bool parallel,int nbCore);
        void perfomOptimPsiHeavyBall(double lipschitzConstant,double gradTol,int nMaxIter,bool parallel,int nbCore);
        void perfomOptimPsiHeavyBallFred(double gradTol,int nMaxIter,bool parallel,int nbCore);
        void perfomOptimPsiHeavyBallGreedy(int nStep,double* step,int nScoring,double* scoring,double gradTol,int nMaxIter,bool parallel,int nbCore);
        void perfomOptimPsiPolyak(int lenAlpha,double* alpha,double gradTol,int nMaxIter,bool parallel,int nbCore);
        void perfomOptimPsiBFGSUnleashDaBeast(double gradTol,int nMaxIter,double w1, double w2, int wMaxIter,int memSize,double leven,bool parallel,int nbCore);
        void perfomOptimPsiBFGSUnleashDaBeastRestart(int lenPsi,double* psi,double gradTol,int nMaxIter,double w1, double w2, int wMaxIter,int memSize,double leven,bool parallel,int nbCore);


        // optimization research
        void perfomOptimPsiBFGSResearch(double gradTol,int nMaxIter,double w1, double w2, int wMaxIter,int memSize,bool parallel,int nbCore);
        void perfomOptimPsiConstantGDResearch(double stepSize,double gradTol,int nMaxIter,bool parallel,int nbCore);
        void perfomOptimPsiGreedyGDResearch(int nStep,double* step,int nScoring,double* scoring,double gradTol,int nMaxIter,bool parallel,int nbCore);
        void perfomOptimPsiStrongWolfGDResearch(double stepSize,double gradTol,int nMaxIter,double w1, double w2, int wMaxIter,bool parallel,int nbCore);
        void perfomOptimPsiBarzilaiBResearch(double gradTol,int nMaxIter,bool parallel,int nbCore);
        void perfomOptimPsiBarzilaiBWolfResearch(double gradTol,double lineSearchRed,double gammaWolf,int nMaxIter,int memSize,bool parallel,int nbCore);
        void perfomOptimPsiBarzilaiBWolf2Research(double gradTol,double lineSearchRed,double gammaWolf,int nMaxIter,int memSize,int maxHiddenIncrease,bool parallel,int nbCore);
        void perfomOptimPsiHeavyBallResearch(double lipschitzConstant,double gradTol,int nMaxIter,bool parallel,int nbCore);
        void perfomOptimPsiHeavyBallFredResearch(double gradTol,int nMaxIter,bool parallel,int nbCore);
        void perfomOptimPsiHeavyBallGreedyResearch(int nStep,double* step,int nScoring,double* scoring,double gradTol,int nMaxIter,bool parallel,int nbCore);
        void perfomOptimPsiPolyakResearch(int lenAlpha,double* alpha,double gradTol,int nMaxIter,bool parallel,int nbCore);
        void perfomOptimPsiBFGSUnleashDaBeastResearch(double gradTol,int nMaxIter,double w1, double w2, int wMaxIter,int memSize,double leven,bool parallel,int nbCore);
        void perfomOptimPsiHessLM(double gradTol,int nMaxIter,double w1, double w2, int wMaxIter,int memSize,double leven,bool parallel,int nbCore);
        

        // brute force searches
        int computeRanking(double cf,double gNorm,int nbHidden,const Eigen::VectorXd& xk, const Eigen::VectorXd& direction, double * tau,int lenTau, double * score);
        int computeRankingCont(double cf,double gNorm,int nbHidden,const Eigen::VectorXd& xk, const Eigen::VectorXd& direction, double * tau,int lenTau, double * score);
        int computeAlphaStar(double cf ,const Eigen::VectorXd& xk, const Eigen::VectorXd& direction, double * alpha,int lenAplha);


        // getters
        void getOptimizedPsi(int nMass,double* mass);
        double getOptimizationtime();
        double getOptimizationCF();
        int getOptimizationCV();
        
        // getters to for research
        std::list<double> getSteps()   {return this->dataConvergence.stp;};
        std::list<double> getNorms()   {return this->dataConvergence.nrm;};
        std::list<double> getTimes()   {return this->dataConvergence.tme;};
        std::list<double> getCostFunc(){return this->dataConvergence.cst;};
        std::list<int> getHidden()     {return this->dataConvergence.hddn;};
        int getnbItt()                 {return this->dataConvergence.nbtt;};
        double getWallTime()           {return this->dataConvergence.wllTme;};

        // GETTER FOR SWIG
        std::list<std::vector<double>> getadjEdge(int);
        std::list<std::vector<double>> getadjPrimalPt(int);
        std::list<std::vector<double>> getCutPosition(int);
        std::list<double> getContributorLine(int);
            
};


#endif