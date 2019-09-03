#include <stdexcept>
#include <iomanip>
#include <limits>
#include "ot3D.h"
#include <sys/time.h>
// #include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/SparseCholesky>
#include <cassert>

////////////////////////////////////////////////////////////////////////////////////////////////
//////////   ____ ___  __  __ ____  _   _ _____ _____    _   _ _____ ___ _     ____  ///////////
//////////  / ___/ _ \|  \/  |  _ \| | | |_   _| ____|  | | | |_   _|_ _| |   / ___| ///////////
////////// | |  | | | | |\/| | |_) | | | | | | |  _|    | | | | | |  | || |   \___ \ ///////////
////////// | |__| |_| | |  | |  __/| |_| | | | | |___   | |_| | | |  | || |___ ___) |///////////
//////////  \____\___/|_|  |_|_|    \___/  |_| |_____|   \___/  |_| |___|_____|____/ ///////////
////////////////////////////////////////////////////////////////////////////////////////////////



void OTproblem::printCGALVERSION(){
    std::cout << "My CGAL library is " <<  CGAL_VERSION_NR << " (1MMmmb1000)" << std::endl; 
    std::cout << std::endl;
    std::cout << "where MM is the major number release, mm is the minor number release, and "
              << "b is the bug fixing number release." << std::endl;
};

void OTproblem::setLaguerrePosition(double* positionX,int nX,double* positionY,int nY,double* positionZ,int nZ){
    if(nX != nY || nX != nZ ){
		throw std::invalid_argument( "Dimension mismatch In Laguerre Position");
    }
    this->Lag.setCoordinates(positionX,positionY,positionZ,nZ);
};
void OTproblem::setLaguerrePsi(double* weight,int nW){ 
    this->nPoints = nW;
    this->Lag.ComputeTesselation3D(weight,nW);
};
void OTproblem::setPolyline(double* positionX,int nX,double* positionY,int nY,double* positionZ,int nZ,double* weight,int nW){
    if(nX != nY || nX != nZ || nX != nW+1){
        throw std::invalid_argument( "Dimension of input arguments of Polyline mismatch" );
    }
    this->pol=Polyline(positionX,positionY,positionZ,weight,nX);   
};

void OTproblem::setPolylineNew(double* positionX,int nX,double* positionY,int nY,double* positionZ,int nZ,double* weight,int nW){
    if(nX != nY || nX != nZ || nX != 2*nW){
        throw std::invalid_argument( "Dimension of input arguments of Polyline mismatch" );
    }
    this->pol=Polyline(positionX,positionY,positionZ,weight,nW,true);   
};

std::list<std::vector<double>> OTproblem::getadjEdge(int i){
    // for( auto l : this->pol.lines){
    //     std::cout<<l.startPoint.x()<<","<<l.startPoint.y()<<std::endl;
    //     std::cout<<l.endPoint.x()<<","<<l.endPoint.y()<<std::endl;
    // }
    std::vector<double> tmp(6);
    std::list<std::vector<double>> ret;
    Rt::Finite_vertices_iterator actualV=this->Lag.T.finite_vertices_begin();
    while ((actualV!=this->Lag.T.finite_vertices_end())&&(actualV->info()!=i)){
        actualV++;
    }
    if (actualV->info()!=i){
        // throw std::invalid_argument( "No such point");
        std::cout<<i<<"? => No such point"<<std::endl;
        return ret;
    }
    
    if (this->Lag.T.dimension()==3){
        std::list<Rt::Facet> facets;
        this->Lag.T.incident_facets(actualV, std::back_inserter(facets));
        K::Segment_3 s;
        K::Ray_3 r;
        CGAL::Object o;
        for (const auto& facet :facets){
            if (!this->Lag.T.is_infinite(facet)){
                o=this->Lag.T.dual(facet);
                if (CGAL::assign(r,o)){
                    s=K::Segment_3(r.source(),r.point(10));
                    tmp[0]=s.source().x();
                    tmp[1]=s.source().y();
                    tmp[2]=s.source().z();
                    tmp[3]=s.target().x();
                    tmp[4]=s.target().y();
                    tmp[5]=s.target().z();               
                }
                if (CGAL::assign(s,o)){
                    tmp[0]=s.source().x();
                    tmp[1]=s.source().y();
                    tmp[2]=s.source().z();
                    tmp[3]=s.target().x();
                    tmp[4]=s.target().y();
                    tmp[5]=s.target().z();
                    ret.push_back(tmp);
                }
            }
        }
    }
    if (this->Lag.T.dimension()==2){
        Point dp0;
        Point dp1;
        std::list<Rt::Facet> facets;
        this->Lag.T.incident_facets(actualV, std::back_inserter(facets));
        int k=0;
        for (const auto& face1 :facets){
            int j=0;
            for (const auto& face2 :facets){
                if ((face1.first->has_neighbor(face2.first))&&(!this->Lag.T.is_infinite(face1))&&(!this->Lag.T.is_infinite(face2))  ) {
                  if ((CGAL::assign(dp0,this->Lag.T.dual(face1)))&&(CGAL::assign(dp1,this->Lag.T.dual(face2)))){
                    tmp[0]=dp0.x();
                    tmp[1]=dp0.y();
                    tmp[2]=dp0.z();
                    tmp[3]=dp1.x();
                    tmp[4]=dp1.y();
                    tmp[5]=dp1.z();
                    ret.push_back(tmp);                    
                  }
                }
                j++;
            }
            k++;
        }
    }
    return ret;
}

std::list<std::vector<double>> OTproblem::getCutPosition( int i ){
    std::list<std::vector<double>> ret;
    Point Pt;
    for (auto& l :pol.lines){
        for (unsigned k=0;k<l.intersectionLocation.size();k++){
            if (l.intersectionLocation[k]==i){
                Pt=l.getPointFromTime(l.intersectionTime[k]);
                ret.push_back({Pt.x(),Pt.y(),Pt.z()});
                Pt=l.getPointFromTime(l.intersectionTime[k+1]);
                ret.push_back({Pt.x(),Pt.y(),Pt.z()});
            }
        }
    }
    return ret;
}



std::list<std::vector<double>> OTproblem::getadjPrimalPt(int i){
    std::list<std::vector<double>> ret;
    std::vector<double> tmp(3);
    for (auto& j : this->Lag.primalPtList[i].adjPrimalPoints){
        tmp[0]=this->Lag.primalPtList[j].x();
        tmp[1]=this->Lag.primalPtList[j].y();
        tmp[2]=this->Lag.primalPtList[j].z();
    }
    ret.push_back(tmp);
    return ret;       
}

void OTproblem::computeAll(int nMass,double* mass,int nCost,double* cost,int nBar,double* bar){
    this->Lag.Intersect(this->pol);
    int i = 0;
    int n=this->Lag.primalPtList.size();
    if (nMass != n  || nCost != n || nBar != 3*n){
        throw std::invalid_argument("Dimensions mismatch in ComputeAll");
    }
    for( auto & primalPt : this->Lag.primalPtList){
        mass[i] = primalPt.mass;
        cost[i] = primalPt.cost;
        bar[3*i] = primalPt.barycenter[0];
        bar[3*i+1] = primalPt.barycenter[1];
        bar[3*i+2] = primalPt.barycenter[2];
        i++;
    }
}

matCOO OTproblem::computeHessian(){
    return this->Lag.Intersex(this->pol);
}

void OTproblem::setMassDiracs(int nMass,double* mass){
    if( nMass != this->nPoints){
        throw std::invalid_argument("Dimensions mismatch in setMassDiracs");
    }
    this-> massDiracs = Eigen::VectorXd(this->nPoints);
    
    for(int i = 0; i < nPoints; i++ ){
        this->massDiracs[i] = mass[i];
    }
}
void OTproblem::computeDeriv(int nMass,double* mass,int nCost,double* cost){
    this->Lag.ComputePolyLineDer(this->pol);
    int i = 0;
    int n=this->pol.lines.size();
    if (nMass != 6*(n)  || nCost != 6*(n)){
        throw std::invalid_argument("Dimensions mismatch in Compute Deriv");
    }   
    // std::abort();
    for (int k=0;k<3*(n+1);k++){
        mass[k]=0.;
        cost[k]=0.;
    }  
    for( const auto & l : this->pol.lines){
        for (int k=0;k<3;k++){
            mass[k+6*i]      = l.massDerivStart[k];
            mass[k+6*i+3]    = l.massDerivEnd[k];
            cost[k+6*i]      = l.costDerivStart[k];
            cost[k+6*i+3]    = l.costDerivEnd[k];
        }
        i++;
    }
}

void OTproblem::computeDerivTotal(int nMass,double* mass,int nCost,double* cost,int nRho,double* rho){
    this->Lag.ComputePolyLineDer(this->pol);
    int i = 0;
    int n=this->pol.lines.size();
    if (nMass != 6*(n)  || nCost != 6*(n)){
        throw std::invalid_argument("Dimensions mismatch in Compute Deriv");
    }   
    // std::abort();
    for (int k=0;k<3*(n+1);k++){
        mass[k]=0.;
        cost[k]=0.;
    }  
    for( const auto & l : this->pol.lines){
        for (int k=0;k<3;k++){
            mass[k+6*i]      = l.massDerivStart[k];
            mass[k+6*i+3]    = l.massDerivEnd[k];
            cost[k+6*i]      = l.costDerivStart[k];
            cost[k+6*i+3]    = l.costDerivEnd[k];
            rho[k+6*i]       = l.rhoDerivStart[k];
            rho[k+6*i+3]     = l.rhoDerivEnd[k];
        }
        i++;
    }
}

void OTproblem::computePolyInfo(int nMass,double* mass,int nCost,double* cost){
    int i = 0;
    int n=this->pol.lines.size();
    if (nMass != n  || nCost != n){
        throw std::invalid_argument("Dimensions mismatch in computeInfo ");
    }
    for( const auto & l : this->pol.lines){
        mass[i]=l.massSeen;
        cost[i]=l.cost;
        i++;
    }
}


double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
derivInfo OTproblem::oracle(Eigen::VectorXd X){
    derivInfo ret;
    this->setLaguerrePsi(X.data(),this->nPoints);
    if(this->parallelComputation){
        this->Lag.IntersectParallel(this->pol,this->nbCore);
    }
    else{
        this->Lag.Intersect(this->pol);
        // this->Lag.Intersect2(this->pol);
    }
    ret.gradient = Eigen::VectorXd(this->nPoints);
    ret.nbHidden = 0;
    Eigen::VectorXd mass = Eigen::VectorXd(this->nPoints);
    Eigen::VectorXd cost = Eigen::VectorXd(this->nPoints);
    int i = 0;
    for( auto & primalPt : this->Lag.primalPtList){
        mass[i] = primalPt.mass;
        if( mass[i] < 1e-12 ) ret.nbHidden++;
        cost[i] = primalPt.cost;
        i++;
    }
    ret.gradient = this->massDiracs - mass;
    ret.costFunction = cost.sum() + X.dot(this->massDiracs -mass);
    std::string colorStart;
    if(ret.nbHidden > this->nPoints/2){
        colorStart = this->red;
    }
    else if(ret.nbHidden > this->nPoints/4){
        colorStart = this->mage;
    }
    else if(ret.nbHidden > this->nPoints/10){
        colorStart = this->yellow;
    }
    else if(ret.nbHidden == 0){
        colorStart = this->blue;
    }
    else{
        colorStart = this->cyan;
    }
    ret.ezHidden = colorStart+std::to_string(ret.nbHidden)+this->reset;
    return ret;
}
///////////////////////////////////////////////////////////////////////////////////////////
///////   ___  ____ _____ ___ __  __    ____  ____   ___   ____ _____ ____ ____   /////////
///////  / _ \|  _ \_   _|_ _|  \/  |  |  _ \|  _ \ / _ \ / ___| ____/ ___/ ___|  /////////
/////// | | | | |_) || |  | || |\/| |  | |_) | |_) | | | | |   |  _| \___ \___ \  /////////
/////// | |_| |  __/ | |  | || |  | |  |  __/|  _ <| |_| | |___| |___ ___) |__) | /////////
///////  \___/|_|    |_| |___|_|  |_|  |_|   |_| \_\\___/ \____|_____|____/____/  /////////
///////////////////////////////////////////////////////////////////////////////////////////                                                                          


void OTproblem::perfomOptimPsiBFGS(double gradTol,int nMaxIter,double w1, double w2, int wMaxIter,int memSize,bool parallel,int nbCore){
    this->parallelComputation = parallel;
    this->nbCore = nbCore;
    double timeExc = get_wall_time();
    int displaySizeit = (int) std::log10(nMaxIter);
    displaySizeit ++;
    int displaySizehid = (int) std::log10(this->nPoints);
    displaySizehid ++;
    Eigen::VectorXd xStart = Eigen::VectorXd(this->nPoints);
    for(int i = 0; i < nPoints; i++){
        xStart[i] = 0.;
    }
    BFGS bf = BFGS((uint)memSize,(uint)this->nPoints);
    int i = 0;
    Eigen::VectorXd d;
    derivInfo current = this->oracle(xStart);
    derivInfo next;
    Eigen::VectorXd x = xStart;
    double g;
    double gradGdotd;
    double s, sm, sp;
    int iterWolfe;
    while(i< nMaxIter && current.gradient.norm() > gradTol){
        // compute step s with wolfe conditions
        d = bf.findDirection(current.gradient);
        s  = 1.0;
        sm =  0.0;
        sp = std::numeric_limits<double>::infinity();
        next = oracle(x+s*d);
        g  = current.costFunction;
        gradGdotd = current.gradient.dot(d);
        iterWolfe = 0;
        while( ((next.costFunction < g + w1*s*gradGdotd) || (next.gradient.dot(d) > w2*gradGdotd)) && iterWolfe < wMaxIter){
            if(next.costFunction < g + w1*s*gradGdotd){
                sp = s;
                s = (sm+sp)/2.;
            }
            else{
                sm=s;
                if(sp == std::numeric_limits<double>::infinity()){
                    s = 2.*s;
                }
                else{
                    s = (sm+sp)/2.;
                }
            }
            next = oracle(x + s*d);
            iterWolfe++;
        }
        std::string decreasecf   = (next.costFunction > g + w1*s*gradGdotd)? blue+"✓"+reset : red+"✘"+reset;
        std::string decreasegrad = (next.gradient.dot(d) < w2*gradGdotd)   ? blue+"✓"+reset : red+"✘"+reset;
        std::cout.precision(5);
        std::cout<<"i = "<<std::setfill('0')<<std::setw(displaySizeit)<<i<<" #Hidden = "<<current.ezHidden<<std::scientific<<" f = "<<current.costFunction<<" ||∇f|| = "<<current.gradient.norm();
        std::cout.precision(2);
        std::cout<<" s = "<<s<<" ↓cf : "<<decreasecf<<" ↓∇f "<<decreasegrad<<std::endl;
        bf.addDirection(s,d,next.gradient,current.gradient);
        current.costFunction = next.costFunction;
        current.gradient     = next.gradient;
        current.nbHidden     = next.nbHidden;
        current.ezHidden     = next.ezHidden;
        x = x + s*d;
        i++;
    }
    timeExc = get_wall_time() -timeExc;
    std::string cv = (i == nMaxIter)? "✘" : "✓";
    std::cout<<"convergence "<<cv<<" Total time "<<std::fixed<<timeExc<<"sec Time per it "<<std::scientific<<timeExc/((double) i)<<"sec"<<std::endl;
    this->optimizedPsi.phiOpt = x;
    this->optimizedPsi.costFunction = current.costFunction;
    this->optimizedPsi.wallTime = timeExc;
    this->optimizedPsi.nbIter = (i == nMaxIter)? -1 : i;
}


void OTproblem::perfomOptimPsiBFGSUnleashDaBeast(double gradTol,int nMaxIter,double w1, double w2, int wMaxIter,int memSize,double leven,bool parallel,int nbCore){
    this->parallelComputation = parallel;
    this->nbCore = nbCore;
    double timeExc = get_wall_time();
    int displaySizeit = (int) std::log10(nMaxIter);
    displaySizeit ++;
    int displaySizehid = (int) std::log10(this->nPoints);
    displaySizehid ++;
    Eigen::VectorXd xStart = Eigen::VectorXd(this->nPoints);
    for(int i = 0; i < nPoints; i++){
        xStart[i] = 0.;
    }
    int numberWithoutHidden=0;
    BFGS bf = BFGS((uint)memSize,(uint)this->nPoints);
    int i = 0;
    Eigen::VectorXd d;
    derivInfo current = this->oracle(xStart);
    derivInfo next;
    Eigen::VectorXd x = xStart;
    double g;
    double gradGdotd;
    double s, sm, sp;
    int iterWolfe;
    while(i< nMaxIter && current.gradient.norm() > gradTol){
        // compute step s with wolfe conditions
        if ( numberWithoutHidden > 1){
            matCOO coucou = this->Lag.Intersex(this->pol);
            Eigen::SparseMatrix<double> Hess =  coucou.toSparseEigen();
            Eigen::SparseMatrix<double> I(coucou.getDimRow(),coucou.getDimCol());
            I.setIdentity();
            I *= leven;
            try{
                Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> chol(Hess-I);
                d = chol.solve(-current.gradient);
            } catch (...){
                std::cout<<"some shit goes wrong"<<std::endl;
                d = bf.findDirection(current.gradient);    
            }
            

        }
        else{
            d = bf.findDirection(current.gradient);
        }
        
        s  = 1.0;
        sm =  0.0;
        sp = std::numeric_limits<double>::infinity();
        next = oracle(x+s*d);
        g  = current.costFunction;
        gradGdotd = current.gradient.dot(d);
        iterWolfe = 0;
        while( ((next.costFunction < g + w1*s*gradGdotd) || (next.gradient.dot(d) > w2*gradGdotd)) && iterWolfe < wMaxIter){
            if(next.costFunction < g + w1*s*gradGdotd){
                sp = s;
                s = (sm+sp)/2.;
            }
            else{
                sm=s;
                if(sp == std::numeric_limits<double>::infinity()){
                    s = 2.*s;
                }
                else{
                    s = (sm+sp)/2.;
                }
            }
            next = oracle(x + s*d);
            iterWolfe++;
        }
        std::string decreasecf   = (next.costFunction > g + w1*s*gradGdotd)? blue+"✓"+reset : red+"✘"+reset;
        std::string decreasegrad = (next.gradient.dot(d) < w2*gradGdotd)   ? blue+"✓"+reset : red+"✘"+reset;
        std::cout.precision(5);
        std::cout<<"i = "<<std::setfill('0')<<std::setw(displaySizeit)<<i<<" #Hidden = "<<current.ezHidden<<std::scientific<<" f = "<<current.costFunction<<" ||∇f|| = "<<current.gradient.norm();
        std::cout.precision(2);
        std::cout<<" s = "<<s<<" ↓cf : "<<decreasecf<<" ↓∇f "<<decreasegrad<<std::endl;
        bf.addDirection(s,d,next.gradient,current.gradient);
        current.costFunction = next.costFunction;
        current.gradient     = next.gradient;
        current.nbHidden     = next.nbHidden;
        current.ezHidden     = next.ezHidden;
        x = x + s*d;
        i++;
        if (current.nbHidden == 0){
            numberWithoutHidden ++;
        }
        else{
            numberWithoutHidden = 0;
        }
    }
    timeExc = get_wall_time() -timeExc;
    std::string cv = (i == nMaxIter)? "✘" : "✓";
    std::cout<<"convergence "<<cv<<" Total time "<<std::fixed<<timeExc<<"sec Time per it "<<std::scientific<<timeExc/((double) i)<<"sec"<<std::endl;
    this->optimizedPsi.phiOpt = x;
    this->optimizedPsi.costFunction = current.costFunction;
    this->optimizedPsi.wallTime = timeExc;
    this->optimizedPsi.nbIter = (i == nMaxIter)? -1 : i;
}

void OTproblem::perfomOptimPsiBFGSUnleashDaBeastRestart(int lenPsi,double* psi,double gradTol,int nMaxIter,double w1, double w2, int wMaxIter,int memSize,double leven,bool parallel,int nbCore){
    this->parallelComputation = parallel;
    this->nbCore = nbCore;
    double timeExc = get_wall_time();
    int displaySizeit = (int) std::log10(nMaxIter);
    displaySizeit ++;
    int displaySizehid = (int) std::log10(this->nPoints);
    displaySizehid ++;
    Eigen::VectorXd xStart = Eigen::VectorXd(this->nPoints);
    for(int i = 0; i < nPoints; i++){
        xStart[i] = psi[i];
    }
    int numberWithoutHidden=0;
    BFGS bf = BFGS((uint)memSize,(uint)this->nPoints);
    int i = 0;
    Eigen::VectorXd d;
    derivInfo current = this->oracle(xStart);
    derivInfo next;
    Eigen::VectorXd x = xStart;
    double g;
    double gradGdotd;
    double s, sm, sp;
    int iterWolfe;
    while(i< nMaxIter && current.gradient.norm() > gradTol){
        // compute step s with wolfe conditions
        if ( numberWithoutHidden > 1){
            matCOO coucou = this->Lag.Intersex(this->pol);
            Eigen::SparseMatrix<double> Hess =  coucou.toSparseEigen();
            Eigen::SparseMatrix<double> I(coucou.getDimRow(),coucou.getDimCol());
            I.setIdentity();
            I *= leven;
            try{
                Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> chol(Hess-I);
                d = chol.solve(-current.gradient);
            } catch (...){
                std::cout<<"some shit goes wrong"<<std::endl;
                d = bf.findDirection(current.gradient);    
            }
            

        }
        else{
            d = bf.findDirection(current.gradient);
        }
        
        s  = 1.0;
        sm =  0.0;
        sp = std::numeric_limits<double>::infinity();
        next = oracle(x+s*d);
        g  = current.costFunction;
        gradGdotd = current.gradient.dot(d);
        iterWolfe = 0;
        while( ((next.costFunction < g + w1*s*gradGdotd) || (next.gradient.dot(d) > w2*gradGdotd)) && iterWolfe < wMaxIter){
            if(next.costFunction < g + w1*s*gradGdotd){
                sp = s;
                s = (sm+sp)/2.;
            }
            else{
                sm=s;
                if(sp == std::numeric_limits<double>::infinity()){
                    s = 2.*s;
                }
                else{
                    s = (sm+sp)/2.;
                }
            }
            next = oracle(x + s*d);
            iterWolfe++;
        }
        std::string decreasecf   = (next.costFunction > g + w1*s*gradGdotd)? blue+"✓"+reset : red+"✘"+reset;
        std::string decreasegrad = (next.gradient.dot(d) < w2*gradGdotd)   ? blue+"✓"+reset : red+"✘"+reset;
        std::cout.precision(5);
        std::cout<<"i = "<<std::setfill('0')<<std::setw(displaySizeit)<<i<<" #Hidden = "<<current.ezHidden<<std::scientific<<" f = "<<current.costFunction<<" ||∇f|| = "<<current.gradient.norm();
        std::cout.precision(2);
        std::cout<<" s = "<<s<<" ↓cf : "<<decreasecf<<" ↓∇f "<<decreasegrad<<std::endl;
        bf.addDirection(s,d,next.gradient,current.gradient);
        current.costFunction = next.costFunction;
        current.gradient     = next.gradient;
        current.nbHidden     = next.nbHidden;
        current.ezHidden     = next.ezHidden;
        x = x + s*d;
        i++;
        if (current.nbHidden == 0){
            numberWithoutHidden ++;
        }
        else{
            numberWithoutHidden = 0;
        }
    }
    timeExc = get_wall_time() -timeExc;
    std::string cv = (i == nMaxIter)? "✘" : "✓";
    std::cout<<"convergence "<<cv<<" Total time "<<std::fixed<<timeExc<<"sec Time per it "<<std::scientific<<timeExc/((double) i)<<"sec"<<std::endl;
    this->optimizedPsi.phiOpt = x;
    this->optimizedPsi.costFunction = current.costFunction;
    this->optimizedPsi.wallTime = timeExc;
    this->optimizedPsi.nbIter = (i == nMaxIter)? -1 : i;
}
void OTproblem::perfomOptimPsiGreedyGD(int nStep,double* step,int nScoring,double* scoring,double gradTol,int nMaxIter,bool parallel,int nbCore){
    this->parallelComputation = parallel;
    this->nbCore = nbCore;
    double timeExc = get_wall_time();
    // fix the size of the output #itt in order to have a super
    // clean display
    int displaySizeit = (int) std::log10(nMaxIter);
    displaySizeit ++;
    // same but for the number of hidden
    int displaySizehid = (int) std::log10(this->nPoints);
    displaySizehid ++;
    Eigen::VectorXd xStart = Eigen::VectorXd(this->nPoints);
    for(int i = 0; i < nPoints; i++){
        xStart[i] = 0.;
    }
    // ici ajouter le paramettre de step dans l'entete
    gradientDescent gD = gradientDescent(-1.0,gradTol,nMaxIter);
    uint i = 0;
    Eigen::VectorXd d;
    derivInfo current = this->oracle(xStart);
    derivInfo next;
    Eigen::VectorXd x = xStart;
    double g;
    double gradGdotd;
    double gradNorm = current.gradient.norm(); 
    double cf = current.costFunction;
    while(i< gD.iMax() && gradNorm > gD.gradN()){
        // compute step s with wolfe conditions
        d = current.gradient;
        gD.setStep(step[computeRankingCont(cf,gradNorm,current.nbHidden,x,current.gradient,step,nStep,scoring)]);
        next = oracle(x+gD.step()*d);
        
        g  = current.costFunction;
        gradGdotd = current.gradient.dot(d);
        // here we test wolf conditions without insuring it is kinda stupid
        std::string decreasecf   = (next.costFunction > g + 0.0001*gD.step()*gradGdotd)? blue+"✓"+reset : red+"✘"+reset;
        std::string decreasegrad = (next.gradient.dot(d) < 0.99*gradGdotd)   ? blue+"✓"+reset : red+"✘"+reset;
        std::cout.precision(5);
        std::cout<<"i = "<<std::setfill('0')<<std::setw(displaySizeit)<<i<<" #Hidden = "<<current.ezHidden<<std::scientific<<" f = "<<current.costFunction<<" ||∇f|| = "<<current.gradient.norm();
        std::cout.precision(2);
        std::cout<<" s = "<<gD.step()<<" ↓cf : "<<decreasecf<<" ↓∇f "<<decreasegrad<<std::endl;
        
        current.costFunction = next.costFunction;
        current.gradient     = next.gradient;
        current.nbHidden     = next.nbHidden;
        current.ezHidden     = next.ezHidden;
        gradNorm = current.gradient.norm();
        cf = current.costFunction;
        x = x + gD.step()*d;
        i++;
    }
    timeExc = get_wall_time() -timeExc;
    std::string cv = (i == gD.iMax())? "✘" : "✓";
    std::cout<<"convergence "<<cv<<" Total time "<<std::fixed<<timeExc<<"sec Time per it "<<std::scientific<<timeExc/((double) i)<<"sec"<<std::endl;
    this->optimizedPsi.phiOpt = x;
    this->optimizedPsi.costFunction = current.costFunction;
    this->optimizedPsi.wallTime = timeExc;
    this->optimizedPsi.nbIter = (i == gD.iMax())? -1 : i;
}

void OTproblem::perfomOptimPsiConstantGD(double stepSize,double gradTol,int nMaxIter,bool parallel,int nbCore){
    this->parallelComputation = parallel;
    this->nbCore = nbCore;
    double timeExc = get_wall_time();
    // fix the size of the output #itt in order to have a super
    // clean display
    int displaySizeit = (int) std::log10(nMaxIter);
    displaySizeit ++;
    // same but for the number of hidden
    int displaySizehid = (int) std::log10(this->nPoints);
    displaySizehid ++;
    Eigen::VectorXd xStart = Eigen::VectorXd(this->nPoints);
    for(int i = 0; i < nPoints; i++){
        xStart[i] = 0.;
    }
    // ici ajouter le paramettre de step dans l'entete
    gradientDescent gD = gradientDescent(stepSize,gradTol,nMaxIter);
    uint i = 0;
    Eigen::VectorXd d;
    derivInfo current = this->oracle(xStart);
    derivInfo next;
    Eigen::VectorXd x = xStart;
    double g;
    double gradGdotd;
    
    while(i< gD.iMax() && current.gradient.norm() > gD.gradN()){
        // compute step s with wolfe conditions
        d = current.gradient;
        next = oracle(x+gD.step()*d);
        
        g  = current.costFunction;
        gradGdotd = current.gradient.dot(d);
        // here we test wolf conditions without insuring it is kinda stupid
        std::string decreasecf   = (next.costFunction > g + 0.0001*gD.step()*gradGdotd)? blue+"✓"+reset : red+"✘"+reset;
        std::string decreasegrad = (next.gradient.dot(d) < 0.99*gradGdotd)   ? blue+"✓"+reset : red+"✘"+reset;
        std::cout.precision(5);
        std::cout<<"i = "<<std::setfill('0')<<std::setw(displaySizeit)<<i<<" #Hidden = "<<current.ezHidden<<std::scientific<<" f = "<<current.costFunction<<" ||∇f|| = "<<current.gradient.norm();
        std::cout.precision(2);
        std::cout<<" s = "<<gD.step()<<" ↓cf : "<<decreasecf<<" ↓∇f "<<decreasegrad<<std::endl;
        
        current.costFunction = next.costFunction;
        current.gradient     = next.gradient;
        current.nbHidden     = next.nbHidden;
        current.ezHidden     = next.ezHidden;
        x = x + gD.step()*d;
        i++;
    }
    timeExc = get_wall_time() -timeExc;
    std::string cv = (i == gD.iMax())? "✘" : "✓";
    std::cout<<"convergence "<<cv<<" Total time "<<std::fixed<<timeExc<<"sec Time per it "<<std::scientific<<timeExc/((double) i)<<"sec"<<std::endl;
    this->optimizedPsi.phiOpt = x;
    this->optimizedPsi.costFunction = current.costFunction;
    this->optimizedPsi.wallTime = timeExc;
    this->optimizedPsi.nbIter = (i == gD.iMax())? -1 : i;
}



void OTproblem::perfomOptimPsiStrongWolfGD(double stepSize,double gradTol,int nMaxIter,double w1, double w2, int wMaxIter,bool parallel,int nbCore){
    this->parallelComputation = parallel;
    this->nbCore = nbCore;
    double timeExc = get_wall_time();
    int displaySizeit = (int) std::log10(nMaxIter);
    displaySizeit ++;
    int displaySizehid = (int) std::log10(this->nPoints);
    displaySizehid ++;
    Eigen::VectorXd xStart = Eigen::VectorXd(this->nPoints);
    for(int i = 0; i < nPoints; i++){
        xStart[i] = 0.;
    }
    gradientStrongWolf gSW = gradientStrongWolf(stepSize, gradTol,w1,w1,wMaxIter,nMaxIter);
    uint i = 0;
    Eigen::VectorXd d;
    derivInfo current = this->oracle(xStart);
    derivInfo next;
    Eigen::VectorXd x = xStart;
    double g;
    double gradGdotd;
    double s, sm, sp;
    uint iterWolfe;
    while(i< gSW.iMax() && current.gradient.norm() > gSW.gradN()){
        // compute step s with wolfe conditions
        d = current.gradient;
        s  = gSW.step();
        sm =  0.0;
        sp = std::numeric_limits<double>::infinity();
        next = oracle(x+s*d);
        g  = current.costFunction;
        gradGdotd = current.gradient.dot(d);
        iterWolfe = 0;
        while( ((next.costFunction < g + gSW.w1()*s*gradGdotd) || (next.gradient.dot(d) > gSW.w2()*gradGdotd)) && iterWolfe < gSW.innerMax()){
            if(next.costFunction < g + w1*s*gradGdotd){
                sp = s;
                s = (sm+sp)/2.;
            }
            else{
                sm=s;
                if(sp == std::numeric_limits<double>::infinity()){
                    s = 2.*s;
                }
                else{
                    s = (sm+sp)/2.;
                }
            }
            next = oracle(x + s*d);
            iterWolfe++;
        }
        std::string decreasecf   = (next.costFunction > g + gSW.w1()*s*gradGdotd)? blue+"✓"+reset : red+"✘"+reset;
        std::string decreasegrad = (next.gradient.dot(d) < gSW.w2()*gradGdotd)   ? blue+"✓"+reset : red+"✘"+reset;
        std::cout.precision(5);
        std::cout<<"i = "<<std::setfill('0')<<std::setw(displaySizeit)<<i<<" #Hidden = "<<current.ezHidden<<std::scientific<<" f = "<<current.costFunction<<" ||∇f|| = "<<current.gradient.norm();
        std::cout.precision(2);
        std::cout<<" s = "<<s<<" ↓cf : "<<decreasecf<<" ↓∇f "<<decreasegrad<<std::endl;
        current.costFunction = next.costFunction;
        current.gradient     = next.gradient;
        current.nbHidden     = next.nbHidden;
        current.ezHidden     = next.ezHidden;
        x = x + s*d;
        i++;
    }
    timeExc = get_wall_time() -timeExc;
    std::string cv = (i == gSW.iMax())? "✘" : "✓";
    std::cout<<"convergence "<<cv<<" Total time "<<std::fixed<<timeExc<<"sec Time per it "<<std::scientific<<timeExc/((double) i)<<"sec"<<std::endl;
    this->optimizedPsi.phiOpt = x;
    this->optimizedPsi.costFunction = current.costFunction;
    this->optimizedPsi.wallTime = timeExc;
    this->optimizedPsi.nbIter = (i == gSW.iMax())? -1 : i;
}



void OTproblem::perfomOptimPsiBarzilaiB(double gradTol,int nMaxIter,bool parallel,int nbCore){
    // ref : Fletcher, R. (2005). On the barzilai-borwein method. 
    // In Optimization and control with applications (pp. 235-256). Springer, Boston, MA.
    this->parallelComputation = parallel;
    this->nbCore = nbCore;
    double timeExc = get_wall_time();
    // fix the size of the output #itt in order to have a super
    // clean display
    int displaySizeit = (int) std::log10(nMaxIter);
    displaySizeit ++;
    // same but for the number of hidden
    int displaySizehid = (int) std::log10(this->nPoints);
    displaySizehid ++;
    Eigen::VectorXd xStart = Eigen::VectorXd(this->nPoints);
    for(int i = 0; i < nPoints; i++){
        xStart[i] = 0.;
    }
    barzilaiBorwein bb = barzilaiBorwein(gradTol,nMaxIter);
    
    
    uint i = 0;
    Eigen::VectorXd dK,dKm1,gammaK,gammaKm1,yKm1;
    double alphaK;
    derivInfo current = this->oracle(xStart);
    derivInfo next;
    Eigen::VectorXd x = xStart;
    while(i< bb.iMax() && current.gradient.norm() > bb.gradN()){
        
        gammaK = - current.gradient;
        if( i == 0){
          alphaK = 10.;
          gammaKm1 = gammaK;
          dK = 0.1*current.gradient;
          dKm1 = 0.1*current.gradient;
        }
        else{
            yKm1 = gammaK - gammaKm1;
            alphaK = dKm1.dot(yKm1)/dKm1.squaredNorm();
            dK = -gammaK/alphaK;
        }

        std::cout.precision(4);
        std::cout<<"i = "<<std::setfill('0')<<std::setw(displaySizeit)<<i<<" #Hidden = "<<current.ezHidden<<std::scientific<<" f = "<<current.costFunction<<" ||∇f|| = "<<current.gradient.norm()<<" s = "<<1.0/alphaK<<std::endl;
        x = x + dK;
        gammaKm1 = gammaK;
        dKm1 = dK;        
        next = oracle(x);
        current.costFunction = next.costFunction;
        current.gradient     = next.gradient;
        current.nbHidden     = next.nbHidden;
        current.ezHidden     = next.ezHidden;
        
        
        
        i++;
    }
    timeExc = get_wall_time() -timeExc;
    std::string cv = (i == bb.iMax())? "✘" : "✓";
    std::cout<<"convergence "<<cv<<" Total time "<<std::fixed<<timeExc<<"sec Time per it "<<std::scientific<<timeExc/((double) i)<<"sec"<<std::endl;
    this->optimizedPsi.phiOpt = x;
    this->optimizedPsi.costFunction = current.costFunction;
    this->optimizedPsi.wallTime = timeExc;
    this->optimizedPsi.nbIter = (i == bb.iMax())? -1 : i;
}


void OTproblem::perfomOptimPsiBarzilaiBWolf(double gradTol,double lineSearchRed,double gammaWolf,int nMaxIter,int memSize,bool parallel,int nbCore){
    // ref : Fletcher, R. (2005). On the barzilai-borwein method. 
    // In Optimization and control with applications (pp. 235-256). Springer, Boston, MA.
    this->parallelComputation = parallel;
    this->nbCore = nbCore;
    double timeExc = get_wall_time();
    // fix the size of the output #itt in order to have a super
    // clean display
    int displaySizeit = (int) std::log10(nMaxIter);
    displaySizeit ++;
    // same but for the number of hidden
    int displaySizehid = (int) std::log10(this->nPoints);
    displaySizehid ++;
    Eigen::VectorXd xStart = Eigen::VectorXd(this->nPoints);
    for(int i = 0; i < nPoints; i++){
        xStart[i] = 0.;
    }
    barzilaiBorweinLS bBG = barzilaiBorweinLS(gradTol,lineSearchRed,gammaWolf,nMaxIter,memSize);
    
    
    uint i = 0;
    int iWolfe = -1; 
    Eigen::VectorXd dK,dKm1,gammaK,gammaKm1,yKm1;
    double alphaK;
    derivInfo current = this->oracle(xStart);
    derivInfo next;
    Eigen::VectorXd x = xStart;
    double lambda;
    derivInfo tmp;
    while(i< bBG.iMax() && current.gradient.norm() > bBG.gradN()){
        bBG.addCost(current.costFunction);
        gammaK = - current.gradient;
        if( i == 0){
          alphaK = 10.;
          gammaKm1 = gammaK;
          dK = 0.1*current.gradient;
          dKm1 = 0.1*current.gradient;
        }
        else{
            yKm1 = gammaK - gammaKm1;
            alphaK = dKm1.dot(yKm1)/dKm1.squaredNorm();
            dK = -gammaK/alphaK;
        }
        if(i > 1){
            lambda = 1/alphaK;
            tmp = oracle(x - lambda*gammaK);
            iWolfe = 0;
            while(tmp.costFunction < bBG.minCost() + bBG.gamma()*lambda*gammaK.squaredNorm()){
                lambda *= bBG.lSearch();
                tmp = oracle(x - lambda*gammaK);
                iWolfe++;
                if( iWolfe > 20){
                    std::cout<<"hum some shit happens gotta do a small gradient ascent step"<<std::endl;
                    lambda = 1e-3;
                    break;
                }
            }            
            dK = -gammaK*lambda;
        }
        

        std::cout.precision(4);
        std::cout<<"i = "<<std::setfill('0')<<std::setw(displaySizeit)<<i<<" #Hidden = "<<current.ezHidden<<std::scientific<<" f = "<<current.costFunction<<" ||∇f|| = "<<current.gradient.norm()<<" s = "<<lambda<<" itWolfe "<<iWolfe<<std::endl;
        x = x + dK;
        gammaKm1 = gammaK;
        dKm1 = dK;        
        next = oracle(x);
        current.costFunction = next.costFunction;
        current.gradient     = next.gradient;
        current.nbHidden     = next.nbHidden;
        current.ezHidden     = next.ezHidden;
        
        
        
        i++;
    }
    timeExc = get_wall_time() -timeExc;
    std::string cv = (i == bBG.iMax())? "✘" : "✓";
    std::cout<<"convergence "<<cv<<" Total time "<<std::fixed<<timeExc<<"sec Time per it "<<std::scientific<<timeExc/((double) i)<<"sec"<<std::endl;
    this->optimizedPsi.phiOpt = x;
    this->optimizedPsi.costFunction = current.costFunction;
    this->optimizedPsi.wallTime = timeExc;
    this->optimizedPsi.nbIter = (i == bBG.iMax())? -1 : i;
}

void OTproblem::perfomOptimPsiBarzilaiBWolf2(double gradTol,double lineSearchRed,double gammaWolf,int nMaxIter,int memSize,int maxHiddenIncrease,bool parallel,int nbCore){
    // ref : Fletcher, R. (2005). On the barzilai-borwein method. 
    // In Optimization and control with applications (pp. 235-256). Springer, Boston, MA.
    this->parallelComputation = parallel;
    this->nbCore = nbCore;
    double timeExc = get_wall_time();
    // fix the size of the output #itt in order to have a super
    // clean display
    int displaySizeit = (int) std::log10(nMaxIter);
    displaySizeit ++;
    // same but for the number of hidden
    int displaySizehid = (int) std::log10(this->nPoints);
    displaySizehid ++;
    Eigen::VectorXd xStart = Eigen::VectorXd(this->nPoints);
    for(int i = 0; i < nPoints; i++){
        xStart[i] = 0.;
    }
    barzilaiBorweinLS bBG = barzilaiBorweinLS(gradTol,lineSearchRed,gammaWolf,nMaxIter,memSize);
    
    uint i = 0;
    int iWolfe = -1; 
    Eigen::VectorXd dK,dKm1,gammaK,gammaKm1,yKm1;
    double alphaK;
    derivInfo current = this->oracle(xStart);
    derivInfo next;
    Eigen::VectorXd x = xStart;
    double lambda;
    derivInfo tmp;
    while(i< bBG.iMax() && current.gradient.norm() > bBG.gradN()){
        bBG.addCost(current.costFunction);
        gammaK = - current.gradient;
        if( i == 0){
          alphaK = 10.;
          gammaKm1 = gammaK;
          dK = 0.1*current.gradient;
          dKm1 = 0.1*current.gradient;
        }
        else{
            yKm1 = gammaK - gammaKm1;
            alphaK = dKm1.dot(yKm1)/dKm1.squaredNorm();
            dK = -gammaK/alphaK;
        }
        if(i > 1){
            lambda = 1/alphaK;
            tmp = oracle(x - lambda*gammaK);
            iWolfe = 0;
            while((tmp.costFunction < bBG.minCost() + bBG.gamma()*lambda*gammaK.squaredNorm()) || ( current.nbHidden + maxHiddenIncrease < tmp.nbHidden)){
                lambda *= bBG.lSearch();
                tmp = oracle(x - lambda*gammaK);
                iWolfe++;
                if( iWolfe > 20){
                    std::cout<<"hum some shit happens gotta do a small gradient ascent step"<<std::endl;
                    lambda = 1e-3;
                    break;
                }
            }            
            dK = -gammaK*lambda;
        }
        

        std::cout.precision(4);
        std::cout<<"i = "<<std::setfill('0')<<std::setw(displaySizeit)<<i<<" #Hidden = "<<current.ezHidden<<std::scientific<<" f = "<<current.costFunction<<" ||∇f|| = "<<current.gradient.norm()<<" s = "<<lambda<<" itWolfe "<<iWolfe<<std::endl;
        x = x + dK;
        gammaKm1 = gammaK;
        dKm1 = dK;        
        next = oracle(x);
        current.costFunction = next.costFunction;
        current.gradient     = next.gradient;
        current.nbHidden     = next.nbHidden;
        current.ezHidden     = next.ezHidden;
        
        
        
        i++;
    }
    timeExc = get_wall_time() -timeExc;
    std::string cv = (i == bBG.iMax())? "✘" : "✓";
    std::cout<<"convergence "<<cv<<" Total time "<<std::fixed<<timeExc<<"sec Time per it "<<std::scientific<<timeExc/((double) i)<<"sec"<<std::endl;
    this->optimizedPsi.phiOpt = x;
    this->optimizedPsi.costFunction = current.costFunction;
    this->optimizedPsi.wallTime = timeExc;
    this->optimizedPsi.nbIter = (i == bBG.iMax())? -1 : i;
}


void OTproblem::perfomOptimPsiHeavyBall(double lipschitzConstant,double gradTol,int nMaxIter,bool parallel,int nbCore){
    this->parallelComputation = parallel;
    this->nbCore = nbCore;
    double timeExc = get_wall_time();
    int displaySizeit = (int) std::log10(nMaxIter);
    displaySizeit ++;
    int displaySizehid = (int) std::log10(this->nPoints);
    displaySizehid ++;
    Eigen::VectorXd xStart = Eigen::VectorXd(this->nPoints);
    for(int i = 0; i < nPoints; i++){
        xStart[i] = 0.;
    }
                    
    heavyBall hBall = heavyBall(lipschitzConstant,gradTol,nMaxIter);
    uint i = 0;
    derivInfo current = this->oracle(xStart);
    derivInfo tmp;
    Eigen::VectorXd yk,xk,xkm1;
    xk = xStart;
    double tau = 1./(hBall.L());
    double cf, gNorm;
    gNorm = current.gradient.norm();
    cf = current.costFunction;
    while(i< hBall.iMax() &&  gNorm > hBall.gradN()){
       if(i ==0){
            xkm1 = xk;
            xk = xk+tau*current.gradient;
        }
        else{
            yk = xk + ((double) i - 1)/((double) i + 2)*(xk-xkm1);
            xkm1 = xk;
            tmp = this->oracle(yk);
            xk = yk + tau*tmp.gradient; 
        }
        current = this->oracle(xk);
        std::cout.precision(5);
        cf = current.costFunction;
        gNorm = current.gradient.norm();
        std::cout<<"i = "<<std::setfill('0')<<std::setw(displaySizeit)<<i<<" #Hidden = "<<current.ezHidden<<std::scientific<<" f = "<<cf<<" ||∇f|| = "<<gNorm<<std::endl;
        i++;
    }
    timeExc = get_wall_time() -timeExc;
    std::string cv = (i == hBall.iMax())? "✘" : "✓";
    std::cout<<"convergence "<<cv<<" Total time "<<std::fixed<<timeExc<<"sec Time per it "<<std::scientific<<timeExc/((double) i)<<"sec"<<std::endl;
    this->optimizedPsi.phiOpt = xk;
    this->optimizedPsi.costFunction = current.costFunction;
    this->optimizedPsi.wallTime = timeExc;
    this->optimizedPsi.nbIter = (i == hBall.iMax())? -1 : i;
}

void OTproblem::perfomOptimPsiHeavyBallFred(double gradTol,int nMaxIter,bool parallel,int nbCore){
    // Nesterov heavy Ball but with the idea on fred to compute empirically an estimate of 
    // the lischitz constant
    this->parallelComputation = parallel;
    this->nbCore = nbCore;
    double timeExc = get_wall_time();
    int displaySizeit = (int) std::log10(nMaxIter);
    displaySizeit ++;
    int displaySizehid = (int) std::log10(this->nPoints);
    displaySizehid ++;
    Eigen::VectorXd xStart = Eigen::VectorXd(this->nPoints);
    for(int i = 0; i < nPoints; i++){
        xStart[i] = 0.;
    }
                    
    heavyBall hBall = heavyBall(-1,gradTol,nMaxIter);
    uint i = 0;
    Eigen::VectorXd d;
    derivInfo current = this->oracle(xStart);
    derivInfo tmp;
    Eigen::VectorXd yk,xk,xkm1,diffGrad,diffx;
    xk = xStart;
    double tau = 0.001;
    double minTau = 100.0;
    while(i< hBall.iMax() && current.gradient.norm() > hBall.gradN()){
        if(i ==0){
            d = current.gradient;
            xkm1 = xk;
            tau = .001;
            xk = xk+tau*d;
        }
        else{
            yk = xk + ((double) i - 1)/((double) i + 2)*(xk-xkm1);
            xkm1 = xk;
            tmp = this->oracle(yk);
            diffGrad = tmp.gradient - current.gradient;
            diffx = yk - xk;
            if ( i > 1){
                tau = diffGrad.norm()/diffx.norm();
                minTau = std::min(tau,minTau);
                tau = minTau;
            }
            else{
                tau = 0.01;
            }            
            xk = yk + tau*tmp.gradient;
        } 
        current = this->oracle(xk);
        std::cout.precision(5);
        std::cout<<"i = "<<std::setfill('0')<<std::setw(displaySizeit)<<i<<" #Hidden = "<<current.ezHidden<<std::scientific<<" f = "<<current.costFunction<<" ||∇f|| = "<<current.gradient.norm()<<" tau = "<<tau<<std::endl;
        i++;
    }
    timeExc = get_wall_time() -timeExc;
    std::string cv = (i == hBall.iMax())? "✘" : "✓";
    std::cout<<"convergence "<<cv<<" Total time "<<std::fixed<<timeExc<<"sec Time per it "<<std::scientific<<timeExc/((double) i)<<"sec"<<std::endl;
    this->optimizedPsi.phiOpt = xk;
    this->optimizedPsi.costFunction = current.costFunction;
    this->optimizedPsi.wallTime = timeExc;
    this->optimizedPsi.nbIter = (i == hBall.iMax())? -1 : i;
}

void OTproblem::perfomOptimPsiHeavyBallGreedy(int nStep,double* step,int nScoring,double* scoring,double gradTol,int nMaxIter,bool parallel,int nbCore){
    this->parallelComputation = parallel;
    this->nbCore = nbCore;
    double timeExc = get_wall_time();
    int displaySizeit = (int) std::log10(nMaxIter);
    displaySizeit ++;
    int displaySizehid = (int) std::log10(this->nPoints);
    displaySizehid ++;
    Eigen::VectorXd xStart = Eigen::VectorXd(this->nPoints);
    for(int i = 0; i < nPoints; i++){
        xStart[i] = 0.;
    }
                    
    heavyBall hBall = heavyBall(-1.0,gradTol,nMaxIter);
    uint i = 0;
    derivInfo current = this->oracle(xStart);
    derivInfo tmp;
    Eigen::VectorXd yk,xk,xkm1;
    xk = xStart;
    double cf, gNorm;
    gNorm = current.gradient.norm();
    cf = current.costFunction;
    int bestTau;
    while(i< hBall.iMax() &&  gNorm > hBall.gradN()){
       if(i ==0){
            xkm1 = xk;
            bestTau = computeRanking(cf,gNorm,current.nbHidden,xk,current.gradient,step,nStep,scoring);
            xk = xk+step[bestTau]*current.gradient;
        }
        else{
            yk = xk + ((double) i - 1)/((double) i + 2)*(xk-xkm1);
            xkm1 = xk;
            tmp = this->oracle(yk);
            bestTau = computeRanking(current.costFunction,current.gradient.norm(),current.nbHidden,yk,tmp.gradient,step,nStep,scoring);
            xk = yk + step[bestTau]*tmp.gradient; 
        }
        current = this->oracle(xk);
        std::cout.precision(5);
        cf = current.costFunction;
        gNorm = current.gradient.norm();
        std::cout<<"i = "<<std::setfill('0')<<std::setw(displaySizeit)<<i<<" #Hidden = "<<current.ezHidden<<std::scientific<<" f = "<<cf<<" ||∇f|| = "<<gNorm<<" tau = "<<step[bestTau]<<std::endl;
        i++;
    }
    timeExc = get_wall_time() -timeExc;
    std::string cv = (i == hBall.iMax())? "✘" : "✓";
    std::cout<<"convergence "<<cv<<" Total time "<<std::fixed<<timeExc<<"sec Time per it "<<std::scientific<<timeExc/((double) i)<<"sec"<<std::endl;
    this->optimizedPsi.phiOpt = xk;
    this->optimizedPsi.costFunction = current.costFunction;
    this->optimizedPsi.wallTime = timeExc;
    this->optimizedPsi.nbIter = (i == hBall.iMax())? -1 : i;
}
void OTproblem::perfomOptimPsiPolyak(int lenAlpha,double* alpha,double gradTol,int nMaxIter,bool parallel,int nbCore){
    // introduction to optimisation Борис Поляк page 72 eqn(24)
    this->parallelComputation = parallel;
    this->nbCore = nbCore;
    double timeExc = get_wall_time();
    int displaySizeit = (int) std::log10(nMaxIter);
    displaySizeit ++;
    int displaySizehid = (int) std::log10(this->nPoints);
    displaySizehid ++;
    Eigen::VectorXd xStart = Eigen::VectorXd(this->nPoints);
    for(int i = 0; i < nPoints; i++){
        xStart[i] = 0.;
    }
                    
    heavyBall hBall = heavyBall(-1.0,gradTol,nMaxIter);
    uint i = 0;
    derivInfo current = this->oracle(xStart);
    derivInfo tmp;
    Eigen::VectorXd xk,pk,pkm1,rk,rkm1;
    xk = xStart;
    double cf, gNorm,beta;
    gNorm = current.gradient.norm();
    cf = current.costFunction;
    int alphaStar;
    while(i< hBall.iMax() &&  gNorm > hBall.gradN()){
       if(i ==0){
            beta = 0;
            rk = current.gradient;
            pk = rk;
            alphaStar = computeAlphaStar(cf,xk,pk,alpha,lenAlpha);
            xk = xk+alpha[alphaStar]*pk;
        }
        else{
            rk = current.gradient;
            beta = rk.dot(rk - rkm1)/rkm1.squaredNorm();
            pk = rk + beta*pkm1;
            alphaStar = computeAlphaStar(cf,xk,pk,alpha,lenAlpha);
            xk = xk+alpha[alphaStar]*pk;
        }
        pkm1 = pk; rkm1 = rk;
        current = this->oracle(xk);
        std::cout.precision(5);
        cf = current.costFunction;
        gNorm = current.gradient.norm();
        std::cout<<"i = "<<std::setfill('0')<<std::setw(displaySizeit)<<i<<" #Hidden = "<<current.ezHidden<<std::scientific<<" f = "<<cf<<" ||∇f|| = "<<gNorm<<" tau = "<<alpha[alphaStar]<<std::endl;
        i++;
    }
    timeExc = get_wall_time() -timeExc;
    std::string cv = (i == hBall.iMax())? "✘" : "✓";
    std::cout<<"convergence "<<cv<<" Total time "<<std::fixed<<timeExc<<"sec Time per it "<<std::scientific<<timeExc/((double) i)<<"sec"<<std::endl;
    this->optimizedPsi.phiOpt = xk;
    this->optimizedPsi.costFunction = current.costFunction;
    this->optimizedPsi.wallTime = timeExc;
    this->optimizedPsi.nbIter = (i == hBall.iMax())? -1 : i;
}


//################################################################################################
//###################### ____  _____ ____  _____    _    ____   ____ _   _   #####################
//######################|  _ \| ____/ ___|| ____|  / \  |  _ \ / ___| | | |  #####################
//######################| |_) |  _| \___ \|  _|   / _ \ | |_) | |   | |_| |  #####################
//######################|  _ <| |___ ___) | |___ / ___ \|  _ <| |___|  _  |  #####################
//######################|_| \_\_____|____/|_____/_/   \_\_| \_\\____|_| |_|  #####################
//######################                                                     #####################
//################################################################################################


void OTproblem::perfomOptimPsiBFGSResearch(double gradTol,int nMaxIter,double w1, double w2, int wMaxIter,int memSize,bool parallel,int nbCore){
    this->parallelComputation = parallel;
    this->nbCore = nbCore;
    double timeExc = get_wall_time();
    int displaySizeit = (int) std::log10(nMaxIter);
    displaySizeit ++;
    int displaySizehid = (int) std::log10(this->nPoints);
    displaySizehid ++;
    Eigen::VectorXd xStart = Eigen::VectorXd(this->nPoints);
    for(int i = 0; i < nPoints; i++){
        xStart[i] = 0.;
    }
    BFGS bf = BFGS((uint)memSize,(uint)this->nPoints);
    int i = 0;
    Eigen::VectorXd d;
    derivInfo current = this->oracle(xStart);
    derivInfo next;
    Eigen::VectorXd x = xStart;
    double g;
    double gradGdotd;
    double s, sm, sp;
    int iterWolfe;
    std::list<double> stploc;
    std::list<double> nrmloc;
    std::list<double> cstloc;
    std::list<double> tmeloc;
    std::list<int> hddnloc;

    nrmloc.push_back(current.gradient.norm());
    cstloc.push_back(current.costFunction);
    hddnloc.push_back(current.nbHidden);
    double timer;
    while(i< nMaxIter && current.gradient.norm() > gradTol){
        timer = get_wall_time();
        // compute step s with wolfe conditions
        d = bf.findDirection(current.gradient);
        s  = 1.0;
        sm =  0.0;
        sp = std::numeric_limits<double>::infinity();
        next = oracle(x+s*d);
        g  = current.costFunction;
        gradGdotd = current.gradient.dot(d);
        iterWolfe = 0;
        while( ((next.costFunction < g + w1*s*gradGdotd) || (next.gradient.dot(d) > w2*gradGdotd)) && iterWolfe < wMaxIter){
            if(next.costFunction < g + w1*s*gradGdotd){
                sp = s;
                s = (sm+sp)/2.;
            }
            else{
                sm=s;
                if(sp == std::numeric_limits<double>::infinity()){
                    s = 2.*s;
                }
                else{
                    s = (sm+sp)/2.;
                }
            }
            next = oracle(x + s*d);
            iterWolfe++;
        }
        stploc.push_back(s);
        nrmloc.push_back(next.gradient.norm());
        cstloc.push_back(next.costFunction);
        hddnloc.push_back(next.nbHidden);

        if ( i % 100 == 0){
            std::string decreasecf   = (next.costFunction > g + w1*s*gradGdotd)? blue+"✓"+reset : red+"✘"+reset;
            std::string decreasegrad = (next.gradient.dot(d) < w2*gradGdotd)   ? blue+"✓"+reset : red+"✘"+reset;
            std::cout.precision(5);
            std::cout<<"i = "<<std::setfill('0')<<std::setw(displaySizeit)<<i<<" #Hidden = "<<current.ezHidden<<std::scientific<<" f = "<<current.costFunction<<" ||∇f|| = "<<current.gradient.norm();
            std::cout.precision(2);
            std::cout<<" s = "<<s<<" ↓cf : "<<decreasecf<<" ↓∇f "<<decreasegrad<<std::endl;
        }
        bf.addDirection(s,d,next.gradient,current.gradient);
        current.costFunction = next.costFunction;
        current.gradient     = next.gradient;
        current.nbHidden     = next.nbHidden;
        current.ezHidden     = next.ezHidden;
        x = x + s*d;
        i++;
        tmeloc.push_back(get_wall_time()-timer);
    }
    timeExc = get_wall_time() -timeExc;
    std::string cv = (i == nMaxIter)? "✘" : "✓";
    std::cout<<"convergence "<<cv<<" Total time "<<std::fixed<<timeExc<<"sec Time per it "<<std::scientific<<timeExc/((double) i)<<"sec"<<std::endl;
    
    this->dataConvergence.stp = stploc;
    this->dataConvergence.nrm = nrmloc;
    this->dataConvergence.cst = cstloc;
    this->dataConvergence.tme = tmeloc;
    this->dataConvergence.hddn = hddnloc;
    this->dataConvergence.nbtt = (i == nMaxIter)? -1 : i;
    this->dataConvergence.wllTme = timeExc;
}


void OTproblem::perfomOptimPsiGreedyGDResearch(int nStep,double* step,int nScoring,double* scoring,double gradTol,int nMaxIter,bool parallel,int nbCore){
    this->parallelComputation = parallel;
    this->nbCore = nbCore;
    std::list<double> stploc;
    std::list<double> nrmloc;
    std::list<double> cstloc;
    std::list<double> tmeloc;
    std::list<int> hddnloc;

    double timeExc = get_wall_time();
    int displaySizeit = (int) std::log10(nMaxIter);
    displaySizeit ++;
    int displaySizehid = (int) std::log10(this->nPoints);
    displaySizehid ++;
    Eigen::VectorXd xStart = Eigen::VectorXd(this->nPoints);
    for(int i = 0; i < nPoints; i++){
        xStart[i] = 0.;
    }
    // ici ajouter le paramettre de step dans l'entete
    gradientDescent gD = gradientDescent(-1.0,gradTol,nMaxIter);
    uint i = 0;
    Eigen::VectorXd d;
    derivInfo current = this->oracle(xStart);
    derivInfo next;
    Eigen::VectorXd x = xStart;
    double g;
    double gradGdotd;
    double gradNorm = current.gradient.norm(); 
    double cf = current.costFunction;
    nrmloc.push_back(gradNorm);
    cstloc.push_back(cf);
    hddnloc.push_back(current.nbHidden);
    double timer;
    while(i< gD.iMax() && gradNorm > gD.gradN()){
        timer = get_wall_time();
        // compute step s with wolfe conditions
        d = current.gradient;
        gD.setStep(step[computeRankingCont(cf,gradNorm,current.nbHidden,x,current.gradient,step,nStep,scoring)]);
        next = oracle(x+gD.step()*d);
        
        if( i % 100 == 0){
            g  = current.costFunction;
            gradGdotd = current.gradient.dot(d);
            // here we test wolf conditions without insuring it is kinda stupid
            std::string decreasecf   = (next.costFunction > g + 0.0001*gD.step()*gradGdotd)? blue+"✓"+reset : red+"✘"+reset;
            std::string decreasegrad = (next.gradient.dot(d) < 0.99*gradGdotd)   ? blue+"✓"+reset : red+"✘"+reset;
            std::cout.precision(5);
            std::cout<<"i = "<<std::setfill('0')<<std::setw(displaySizeit)<<i<<" #Hidden = "<<current.ezHidden<<std::scientific<<" f = "<<current.costFunction<<" ||∇f|| = "<<current.gradient.norm();
            std::cout.precision(2);
            std::cout<<" s = "<<gD.step()<<" ↓cf : "<<decreasecf<<" ↓∇f "<<decreasegrad<<std::endl;
        }
        current.costFunction = next.costFunction;
        current.gradient     = next.gradient;
        current.nbHidden     = next.nbHidden;
        current.ezHidden     = next.ezHidden;
        gradNorm = current.gradient.norm();
        cf = current.costFunction;
        stploc.push_back(gD.step());
        nrmloc.push_back(gradNorm);
        cstloc.push_back(cf);
        hddnloc.push_back(current.nbHidden);
        x = x + gD.step()*d;
        i++;
        tmeloc.push_back(get_wall_time()-timer);
    }
    timeExc = get_wall_time() -timeExc;
    std::string cv = (i == gD.iMax())? "✘" : "✓";
    std::cout<<"convergence "<<cv<<" Total time "<<std::fixed<<timeExc<<"sec Time per it "<<std::scientific<<timeExc/((double) i)<<"sec"<<std::endl;
    this->dataConvergence.stp = stploc;
    this->dataConvergence.nrm = nrmloc;
    this->dataConvergence.cst = cstloc;
    this->dataConvergence.tme = tmeloc;
    this->dataConvergence.hddn = hddnloc;
    this->dataConvergence.nbtt = ((int) i == nMaxIter)? -1 : (int) i;
    this->dataConvergence.wllTme = timeExc;

}

void OTproblem::perfomOptimPsiConstantGDResearch(double stepSize,double gradTol,int nMaxIter,bool parallel,int nbCore){
    this->parallelComputation = parallel;
    this->nbCore = nbCore;
    std::list<double> stploc;
    std::list<double> nrmloc;
    std::list<double> cstloc;
    std::list<double> tmeloc;
    std::list<int> hddnloc; 
    double timeExc = get_wall_time();
    // fix the size of the output #itt in order to have a super
    // clean display
    int displaySizeit = (int) std::log10(nMaxIter);
    displaySizeit ++;
    // same but for the number of hidden
    int displaySizehid = (int) std::log10(this->nPoints);
    displaySizehid ++;
    Eigen::VectorXd xStart = Eigen::VectorXd(this->nPoints);
    for(int i = 0; i < nPoints; i++){
        xStart[i] = 0.;
    }
    // ici ajouter le paramettre de step dans l'entete
    gradientDescent gD = gradientDescent(stepSize,gradTol,nMaxIter);
    uint i = 0;
    Eigen::VectorXd d;
    derivInfo current = this->oracle(xStart);
    derivInfo next;
    Eigen::VectorXd x = xStart;
    double g;
    double gradGdotd;
    double cf, gradNorm;
    cf = current.costFunction;
    gradNorm = current.gradient.norm();
    nrmloc.push_back(gradNorm);
    cstloc.push_back(cf);
    hddnloc.push_back(current.nbHidden);
    double timer;
    while(i< gD.iMax() && gradNorm > gD.gradN()){
        timer = get_wall_time();
        // compute step s with wolfe conditions
        d = current.gradient;
        next = oracle(x+gD.step()*d);
        
        g  = current.costFunction;
        gradGdotd = current.gradient.dot(d);
        // here we test wolf conditions without insuring it is kinda stupid
        if( i % 100 == 0){
            std::string decreasecf   = (next.costFunction > g + 0.0001*gD.step()*gradGdotd)? blue+"✓"+reset : red+"✘"+reset;
            std::string decreasegrad = (next.gradient.dot(d) < 0.99*gradGdotd)   ? blue+"✓"+reset : red+"✘"+reset;
            std::cout.precision(5);
            std::cout<<"i = "<<std::setfill('0')<<std::setw(displaySizeit)<<i<<" #Hidden = "<<current.ezHidden<<std::scientific<<" f = "<<current.costFunction<<" ||∇f|| = "<<current.gradient.norm();
            std::cout.precision(2);
            std::cout<<" s = "<<gD.step()<<" ↓cf : "<<decreasecf<<" ↓∇f "<<decreasegrad<<std::endl;
        }


        current.costFunction = next.costFunction;
        current.gradient     = next.gradient;
        current.nbHidden     = next.nbHidden;
        current.ezHidden     = next.ezHidden;

        gradNorm = current.gradient.norm();
        cf = current.costFunction;
        stploc.push_back(gD.step());
        nrmloc.push_back(gradNorm);
        cstloc.push_back(cf);
        hddnloc.push_back(current.nbHidden);

        x = x + gD.step()*d;
        i++;
        tmeloc.push_back(get_wall_time()-timer);
    }
    timeExc = get_wall_time() -timeExc;
    std::string cv = (i == gD.iMax())? "✘" : "✓";
    std::cout<<"convergence "<<cv<<" Total time "<<std::fixed<<timeExc<<"sec Time per it "<<std::scientific<<timeExc/((double) i)<<"sec"<<std::endl;
    this->dataConvergence.stp = stploc;
    this->dataConvergence.nrm = nrmloc;
    this->dataConvergence.cst = cstloc;
    this->dataConvergence.tme = tmeloc;
    this->dataConvergence.hddn = hddnloc;
    this->dataConvergence.nbtt = ((int) i == nMaxIter)? -1 : (int) i;
    this->dataConvergence.wllTme = timeExc;
}

void OTproblem::perfomOptimPsiStrongWolfGDResearch(double stepSize,double gradTol,int nMaxIter,double w1, double w2, int wMaxIter,bool parallel,int nbCore){
    this->parallelComputation = parallel;
    this->nbCore = nbCore;
    std::list<double> stploc;
    std::list<double> nrmloc;
    std::list<double> cstloc;
    std::list<double> tmeloc;
    std::list<int> hddnloc;
    double timeExc = get_wall_time();
    int displaySizeit = (int) std::log10(nMaxIter);
    displaySizeit ++;
    int displaySizehid = (int) std::log10(this->nPoints);
    displaySizehid ++;
    Eigen::VectorXd xStart = Eigen::VectorXd(this->nPoints);
    for(int i = 0; i < nPoints; i++){
        xStart[i] = 0.;
    }
    gradientStrongWolf gSW = gradientStrongWolf(stepSize, gradTol,w1,w1,wMaxIter,nMaxIter);
    uint i = 0;
    Eigen::VectorXd d;
    derivInfo current = this->oracle(xStart);
    derivInfo next;
    Eigen::VectorXd x = xStart;
    double g;
    double gradGdotd;
    double s, sm, sp;
    uint iterWolfe;
    double cf, gradNorm;
    cf = current.costFunction;
    gradNorm = current.gradient.norm();
    nrmloc.push_back(gradNorm);
    cstloc.push_back(cf);
    hddnloc.push_back(current.nbHidden);
    double timer;
    
    while(i< gSW.iMax() && gradNorm > gSW.gradN()){
        timer = get_wall_time();
        // compute step s with wolfe conditions
        d = current.gradient;
        s  = gSW.step();
        sm =  0.0;
        sp = std::numeric_limits<double>::infinity();
        next = oracle(x+s*d);
        g  = current.costFunction;
        gradGdotd = current.gradient.dot(d);
        iterWolfe = 0;
        while( ((next.costFunction < g + gSW.w1()*s*gradGdotd) || (next.gradient.dot(d) > gSW.w2()*gradGdotd)) && iterWolfe < gSW.innerMax()){
            if(next.costFunction < g + w1*s*gradGdotd){
                sp = s;
                s = (sm+sp)/2.;
            }
            else{
                sm=s;
                if(sp == std::numeric_limits<double>::infinity()){
                    s = 2.*s;
                }
                else{
                    s = (sm+sp)/2.;
                }
            }
            next = oracle(x + s*d);
            iterWolfe++;
        }

        if( i % 100 == 0){
            std::string decreasecf   = (next.costFunction > g + gSW.w1()*s*gradGdotd)? blue+"✓"+reset : red+"✘"+reset;
            std::string decreasegrad = (next.gradient.dot(d) < gSW.w2()*gradGdotd)   ? blue+"✓"+reset : red+"✘"+reset;
            std::cout.precision(5);
            std::cout<<"i = "<<std::setfill('0')<<std::setw(displaySizeit)<<i<<" #Hidden = "<<current.ezHidden<<std::scientific<<" f = "<<current.costFunction<<" ||∇f|| = "<<current.gradient.norm();
            std::cout.precision(2);
            std::cout<<" s = "<<s<<" ↓cf : "<<decreasecf<<" ↓∇f "<<decreasegrad<<std::endl;
        }
        current.costFunction = next.costFunction;
        current.gradient     = next.gradient;
        current.nbHidden     = next.nbHidden;
        current.ezHidden     = next.ezHidden;
        
        gradNorm = current.gradient.norm();
        cf = current.costFunction;
        stploc.push_back(s);
        nrmloc.push_back(gradNorm);
        cstloc.push_back(cf);
        hddnloc.push_back(current.nbHidden);
        x = x + s*d;
        i++;
        tmeloc.push_back(get_wall_time()-timer);
    }
    timeExc = get_wall_time() -timeExc;
    std::string cv = (i == gSW.iMax())? "✘" : "✓";
    std::cout<<"convergence "<<cv<<" Total time "<<std::fixed<<timeExc<<"sec Time per it "<<std::scientific<<timeExc/((double) i)<<"sec"<<std::endl;
    this->dataConvergence.stp = stploc;
    this->dataConvergence.nrm = nrmloc;
    this->dataConvergence.cst = cstloc;
    this->dataConvergence.tme = tmeloc;
    this->dataConvergence.hddn = hddnloc;
    this->dataConvergence.nbtt = ((int) i == nMaxIter)? -1 : (int) i;
    this->dataConvergence.wllTme = timeExc;
}



void OTproblem::perfomOptimPsiBarzilaiBResearch(double gradTol,int nMaxIter,bool parallel,int nbCore){
    // ref : Fletcher, R. (2005). On the barzilai-borwein method. 
    // In Optimization and control with applications (pp. 235-256). Springer, Boston, MA.
    this->parallelComputation = parallel;
    this->nbCore = nbCore;
    std::list<double> stploc;
    std::list<double> nrmloc;
    std::list<double> cstloc;
    std::list<double> tmeloc;
    std::list<int> hddnloc;
    double timeExc = get_wall_time();
    // fix the size of the output #itt in order to have a super
    // clean display
    int displaySizeit = (int) std::log10(nMaxIter);
    displaySizeit ++;
    // same but for the number of hidden
    int displaySizehid = (int) std::log10(this->nPoints);
    displaySizehid ++;
    Eigen::VectorXd xStart = Eigen::VectorXd(this->nPoints);
    for(int i = 0; i < nPoints; i++){
        xStart[i] = 0.;
    }
    barzilaiBorwein bb = barzilaiBorwein(gradTol,nMaxIter);
    
    
    uint i = 0;
    Eigen::VectorXd dK,dKm1,gammaK,gammaKm1,yKm1;
    double alphaK;
    derivInfo current = this->oracle(xStart);
    derivInfo next;
    Eigen::VectorXd x = xStart;
    
    double cf, gradNorm;
    cf = current.costFunction;
    gradNorm = current.gradient.norm();
    nrmloc.push_back(gradNorm);
    cstloc.push_back(cf);
    hddnloc.push_back(current.nbHidden);
    double timer;

    while(i< bb.iMax() && gradNorm > bb.gradN()){
        timer = get_wall_time();
        gammaK = - current.gradient;
        if( i == 0){
          alphaK = 10.;
          gammaKm1 = gammaK;
          dK = 0.1*current.gradient;
          dKm1 = 0.1*current.gradient;
        }
        else{
            yKm1 = gammaK - gammaKm1;
            alphaK = dKm1.dot(yKm1)/dKm1.squaredNorm();
            dK = -gammaK/alphaK;
        }

        if( i % 100 == 0){
            std::cout.precision(4);
            std::cout<<"i = "<<std::setfill('0')<<std::setw(displaySizeit)<<i<<" #Hidden = "<<current.ezHidden<<std::scientific<<" f = "<<current.costFunction<<" ||∇f|| = "<<current.gradient.norm()<<" s = "<<1.0/alphaK<<std::endl;
        }

        x = x + dK;
        gammaKm1 = gammaK;
        dKm1 = dK;        
        next = oracle(x);
        current.costFunction = next.costFunction;
        current.gradient     = next.gradient;
        current.nbHidden     = next.nbHidden;
        current.ezHidden     = next.ezHidden;

        gradNorm = current.gradient.norm();
        cf = current.costFunction;
        stploc.push_back(1./alphaK);
        nrmloc.push_back(gradNorm);
        cstloc.push_back(cf);
        hddnloc.push_back(current.nbHidden);
        
        i++;
        tmeloc.push_back(get_wall_time()-timer);
    }
    timeExc = get_wall_time() -timeExc;
    std::string cv = (i == bb.iMax())? "✘" : "✓";
    std::cout<<"convergence "<<cv<<" Total time "<<std::fixed<<timeExc<<"sec Time per it "<<std::scientific<<timeExc/((double) i)<<"sec"<<std::endl;
    this->dataConvergence.stp = stploc;
    this->dataConvergence.nrm = nrmloc;
    this->dataConvergence.cst = cstloc;
    this->dataConvergence.tme = tmeloc;
    this->dataConvergence.hddn = hddnloc;
    this->dataConvergence.nbtt = ((int) i == nMaxIter)? -1 : (int) i;
    this->dataConvergence.wllTme = timeExc;
}


void OTproblem::perfomOptimPsiBarzilaiBWolfResearch(double gradTol,double lineSearchRed,double gammaWolf,int nMaxIter,int memSize,bool parallel,int nbCore){
    // ref : Fletcher, R. (2005). On the barzilai-borwein method. 
    // In Optimization and control with applications (pp. 235-256). Springer, Boston, MA.
    this->parallelComputation = parallel;
    this->nbCore = nbCore;
    std::list<double> stploc;
    std::list<double> nrmloc;
    std::list<double> cstloc;
    std::list<double> tmeloc;
    std::list<int> hddnloc;
    double timeExc = get_wall_time();
    // fix the size of the output #itt in order to have a super
    // clean display
    int displaySizeit = (int) std::log10(nMaxIter);
    displaySizeit ++;
    // same but for the number of hidden
    int displaySizehid = (int) std::log10(this->nPoints);
    displaySizehid ++;
    Eigen::VectorXd xStart = Eigen::VectorXd(this->nPoints);
    for(int i = 0; i < nPoints; i++){
        xStart[i] = 0.;
    }
    barzilaiBorweinLS bBG = barzilaiBorweinLS(gradTol,lineSearchRed,gammaWolf,nMaxIter,memSize);
    
    
    uint i = 0;
    int iWolfe = -1; 
    Eigen::VectorXd dK,dKm1,gammaK,gammaKm1,yKm1;
    double alphaK;
    derivInfo current = this->oracle(xStart);
    derivInfo next;
    Eigen::VectorXd x = xStart;
    double lambda;
    derivInfo tmp;
    double cf, gradNorm;
    cf = current.costFunction;
    gradNorm = current.gradient.norm();
    nrmloc.push_back(gradNorm);
    cstloc.push_back(cf);
    hddnloc.push_back(current.nbHidden);
    double timer;

    while(i< bBG.iMax() && gradNorm > bBG.gradN()){
        timer = get_wall_time();
        bBG.addCost(current.costFunction);
        gammaK = - current.gradient;
        if( i == 0){
          alphaK = 10.;
          gammaKm1 = gammaK;
          dK = 0.1*current.gradient;
          dKm1 = 0.1*current.gradient;
          lambda = 0.1;
        }
        else{
            yKm1 = gammaK - gammaKm1;
            alphaK = dKm1.dot(yKm1)/dKm1.squaredNorm();
            lambda = 1/alphaK;
            dK = -gammaK/alphaK;
        }
        if(i > 1){
            lambda = 1/alphaK;
            tmp = oracle(x - lambda*gammaK);
            iWolfe = 0;
            while(tmp.costFunction < bBG.minCost() + bBG.gamma()*lambda*gammaK.squaredNorm()){
                lambda *= bBG.lSearch();
                tmp = oracle(x - lambda*gammaK);
                iWolfe++;
                if( iWolfe > 20){
                    std::cout<<"hum some shit happens gotta do a small gradient ascent step"<<std::endl;
                    lambda = 1e-3;
                    break;
                }
            }            
            dK = -gammaK*lambda;
        }
        
        if( i % 100 == 0){
            std::cout.precision(4);
            std::cout<<"i = "<<std::setfill('0')<<std::setw(displaySizeit)<<i<<" #Hidden = "<<current.ezHidden<<std::scientific<<" f = "<<current.costFunction<<" ||∇f|| = "<<current.gradient.norm()<<" s = "<<lambda<<" itWolfe "<<iWolfe<<std::endl;
        }

        x = x + dK;
        gammaKm1 = gammaK;
        dKm1 = dK;        
        next = oracle(x);
        current.costFunction = next.costFunction;
        current.gradient     = next.gradient;
        current.nbHidden     = next.nbHidden;
        current.ezHidden     = next.ezHidden;
        
        gradNorm = current.gradient.norm();
        cf = current.costFunction;
        stploc.push_back(lambda);
        nrmloc.push_back(gradNorm);
        cstloc.push_back(cf);
        hddnloc.push_back(current.nbHidden);
        
        i++;
        tmeloc.push_back(get_wall_time()-timer);
    }
    timeExc = get_wall_time() -timeExc;
    std::string cv = (i == bBG.iMax())? "✘" : "✓";
    std::cout<<"convergence "<<cv<<" Total time "<<std::fixed<<timeExc<<"sec Time per it "<<std::scientific<<timeExc/((double) i)<<"sec"<<std::endl;
    this->dataConvergence.stp = stploc;
    this->dataConvergence.nrm = nrmloc;
    this->dataConvergence.cst = cstloc;
    this->dataConvergence.tme = tmeloc;
    this->dataConvergence.hddn = hddnloc;
    this->dataConvergence.nbtt = ((int) i == nMaxIter)? -1 : (int) i;
    this->dataConvergence.wllTme = timeExc;
}


void OTproblem::perfomOptimPsiBarzilaiBWolf2Research(double gradTol,double lineSearchRed,double gammaWolf,int nMaxIter,int memSize,int maxHiddenIncrease,bool parallel,int nbCore){
    // ref : Fletcher, R. (2005). On the barzilai-borwein method. 
    // In Optimization and control with applications (pp. 235-256). Springer, Boston, MA.
    this->parallelComputation = parallel;
    this->nbCore = nbCore;
    std::list<double> stploc;
    std::list<double> nrmloc;
    std::list<double> cstloc;
    std::list<double> tmeloc;
    std::list<int> hddnloc;
    double timeExc = get_wall_time();    
    int displaySizeit = (int) std::log10(nMaxIter);
    displaySizeit ++;

    int displaySizehid = (int) std::log10(this->nPoints);
    displaySizehid ++;
    Eigen::VectorXd xStart = Eigen::VectorXd(this->nPoints);
    for(int i = 0; i < nPoints; i++){
        xStart[i] = 0.;
    }
    barzilaiBorweinLS bBG = barzilaiBorweinLS(gradTol,lineSearchRed,gammaWolf,nMaxIter,memSize);
    
    uint i = 0;
    int iWolfe = -1; 
    Eigen::VectorXd dK,dKm1,gammaK,gammaKm1,yKm1;
    double alphaK;
    derivInfo current = this->oracle(xStart);
    derivInfo next;
    Eigen::VectorXd x = xStart;
    double lambda;
    derivInfo tmp;

    double cf, gradNorm;
    cf = current.costFunction;
    gradNorm = current.gradient.norm();
    nrmloc.push_back(gradNorm);
    cstloc.push_back(cf);
    hddnloc.push_back(current.nbHidden);
    double timer;

    while(i< bBG.iMax() && gradNorm > bBG.gradN()){
        timer = get_wall_time();
        bBG.addCost(current.costFunction);
        gammaK = - current.gradient;
        if( i == 0){
          alphaK = 10.;
          gammaKm1 = gammaK;
          dK = 0.1*current.gradient;
          dKm1 = 0.1*current.gradient;
          lambda = 0.1;
        }
        else{
            yKm1 = gammaK - gammaKm1;
            alphaK = dKm1.dot(yKm1)/dKm1.squaredNorm();
            lambda = 1.0/alphaK;
            dK = -gammaK/alphaK;
        }
        if(i > 1){
            lambda = 1/alphaK;
            tmp = oracle(x - lambda*gammaK);
            iWolfe = 0;
            while((tmp.costFunction < bBG.minCost() + bBG.gamma()*lambda*gammaK.squaredNorm()) || ( current.nbHidden + maxHiddenIncrease < tmp.nbHidden)){
                lambda *= bBG.lSearch();
                tmp = oracle(x - lambda*gammaK);
                iWolfe++;
                if( iWolfe > 20){
                    std::cout<<"hum some shit happens gotta do a small gradient ascent step"<<std::endl;
                    lambda = 1e-3;
                    break;
                }
            }            
            dK = -gammaK*lambda;
        }
        
        if( i % 100 == 0){
            std::cout.precision(4);
            std::cout<<"i = "<<std::setfill('0')<<std::setw(displaySizeit)<<i<<" #Hidden = "<<current.ezHidden<<std::scientific<<" f = "<<current.costFunction<<" ||∇f|| = "<<current.gradient.norm()<<" s = "<<lambda<<" itWolfe "<<iWolfe<<std::endl;
        }
        x = x + dK;

        gammaKm1 = gammaK;
        dKm1 = dK;        
        next = oracle(x);
        
        current.costFunction = next.costFunction;
        current.gradient     = next.gradient;
        current.nbHidden     = next.nbHidden;
        current.ezHidden     = next.ezHidden;
        
        gradNorm = current.gradient.norm();
        cf = current.costFunction;
        stploc.push_back(lambda);
        nrmloc.push_back(gradNorm);
        cstloc.push_back(cf);
        hddnloc.push_back(current.nbHidden);
        
        i++;
        tmeloc.push_back(get_wall_time()-timer);
    }
    timeExc = get_wall_time() -timeExc;
    std::string cv = (i == bBG.iMax())? "✘" : "✓";
    std::cout<<"convergence "<<cv<<" Total time "<<std::fixed<<timeExc<<"sec Time per it "<<std::scientific<<timeExc/((double) i)<<"sec"<<std::endl;
    this->dataConvergence.stp = stploc;
    this->dataConvergence.nrm = nrmloc;
    this->dataConvergence.cst = cstloc;
    this->dataConvergence.tme = tmeloc;
    this->dataConvergence.hddn = hddnloc;
    this->dataConvergence.nbtt = ((int) i == nMaxIter)? -1 : (int) i;
    this->dataConvergence.wllTme = timeExc;
}

void OTproblem::perfomOptimPsiHeavyBallResearch(double lipschitzConstant,double gradTol,int nMaxIter,bool parallel,int nbCore){
    this->parallelComputation = parallel;
    this->nbCore = nbCore;
    std::list<double> stploc;
    std::list<double> nrmloc;
    std::list<double> cstloc;
    std::list<double> tmeloc;
    std::list<int> hddnloc;

    double timeExc = get_wall_time();
    int displaySizeit = (int) std::log10(nMaxIter);
    displaySizeit ++;
    int displaySizehid = (int) std::log10(this->nPoints);
    displaySizehid ++;
    Eigen::VectorXd xStart = Eigen::VectorXd(this->nPoints);
    for(int i = 0; i < nPoints; i++){
        xStart[i] = 0.;
    }
                    
    heavyBall hBall = heavyBall(lipschitzConstant,gradTol,nMaxIter);
    uint i = 0;
    derivInfo current = this->oracle(xStart);
    derivInfo tmp;
    Eigen::VectorXd yk,xk,xkm1;
    xk = xStart;
    double tau = 1./(hBall.L());
    double cf, gNorm;
    gNorm = current.gradient.norm();
    cf = current.costFunction;
    nrmloc.push_back(gNorm);
    cstloc.push_back(cf);
    hddnloc.push_back(current.nbHidden);
    double timer;

    while(i< hBall.iMax() &&  gNorm > hBall.gradN()){
        timer = get_wall_time();
       if(i ==0){
            xkm1 = xk;
            xk = xk+tau*current.gradient;
        }
        else{
            yk = xk + ((double) i - 1)/((double) i + 2)*(xk-xkm1);
            xkm1 = xk;
            tmp = this->oracle(yk);
            xk = yk + tau*tmp.gradient; 
        }
        current = this->oracle(xk);
        std::cout.precision(5);
        cf = current.costFunction;
        gNorm = current.gradient.norm();
        stploc.push_back(tau);
        nrmloc.push_back(gNorm);
        cstloc.push_back(cf);
        hddnloc.push_back(current.nbHidden);
        if( i % 100 == 0){
            std::cout<<"i = "<<std::setfill('0')<<std::setw(displaySizeit)<<i<<" #Hidden = "<<current.ezHidden<<std::scientific<<" f = "<<cf<<" ||∇f|| = "<<gNorm<<std::endl;
        }
        
        i++;
        tmeloc.push_back(get_wall_time()-timer);
    }
    timeExc = get_wall_time() -timeExc;
    std::string cv = (i == hBall.iMax())? "✘" : "✓";
    std::cout<<"convergence "<<cv<<" Total time "<<std::fixed<<timeExc<<"sec Time per it "<<std::scientific<<timeExc/((double) i)<<"sec"<<std::endl;
    this->dataConvergence.stp = stploc;
    this->dataConvergence.nrm = nrmloc;
    this->dataConvergence.cst = cstloc;
    this->dataConvergence.tme = tmeloc;
    this->dataConvergence.hddn = hddnloc;
    this->dataConvergence.nbtt = ((int) i == nMaxIter)? -1 : (int) i;
    this->dataConvergence.wllTme = timeExc;
}

void OTproblem::perfomOptimPsiHeavyBallFredResearch(double gradTol,int nMaxIter,bool parallel,int nbCore){
    // Nesterov heavy Ball but with the idea on fred to compute empirically an estimate of 
    // the lischitz constant
    this->parallelComputation = parallel;
    this->nbCore = nbCore;
    double timeExc = get_wall_time();
    int displaySizeit = (int) std::log10(nMaxIter);
    displaySizeit ++;
    int displaySizehid = (int) std::log10(this->nPoints);
    displaySizehid ++;
    Eigen::VectorXd xStart = Eigen::VectorXd(this->nPoints);
    for(int i = 0; i < nPoints; i++){
        xStart[i] = 0.;
    }
                    
    heavyBall hBall = heavyBall(-1,gradTol,nMaxIter);
    uint i = 0;
    Eigen::VectorXd d;
    derivInfo current = this->oracle(xStart);
    derivInfo tmp;
    Eigen::VectorXd yk,xk,xkm1,diffGrad,diffx;
    xk = xStart;
    double tau = 0.001;
    double minTau = 100.0;
    while(i< hBall.iMax() && current.gradient.norm() > hBall.gradN()){
        if(i ==0){
            d = current.gradient;
            xkm1 = xk;
            tau = .001;
            xk = xk+tau*d;
        }
        else{
            yk = xk + ((double) i - 1)/((double) i + 2)*(xk-xkm1);
            xkm1 = xk;
            tmp = this->oracle(yk);
            diffGrad = tmp.gradient - current.gradient;
            diffx = yk - xk;
            if ( i > 1){
                tau = diffGrad.norm()/diffx.norm();
                minTau = std::min(tau,minTau);
                tau = minTau;
            }
            else{
                tau = 0.01;
            }            
            xk = yk + tau*tmp.gradient;
        } 
        current = this->oracle(xk);
        std::cout.precision(5);
        std::cout<<"i = "<<std::setfill('0')<<std::setw(displaySizeit)<<i<<" #Hidden = "<<current.ezHidden<<std::scientific<<" f = "<<current.costFunction<<" ||∇f|| = "<<current.gradient.norm()<<" tau = "<<tau<<std::endl;
        i++;
    }
    timeExc = get_wall_time() -timeExc;
    std::string cv = (i == hBall.iMax())? "✘" : "✓";
    std::cout<<"convergence "<<cv<<" Total time "<<std::fixed<<timeExc<<"sec Time per it "<<std::scientific<<timeExc/((double) i)<<"sec"<<std::endl;
    this->optimizedPsi.phiOpt = xk;
    this->optimizedPsi.costFunction = current.costFunction;
    this->optimizedPsi.wallTime = timeExc;
    this->optimizedPsi.nbIter = (i == hBall.iMax())? -1 : i;
}

void OTproblem::perfomOptimPsiHeavyBallGreedyResearch(int nStep,double* step,int nScoring,double* scoring,double gradTol,int nMaxIter,bool parallel,int nbCore){
    this->parallelComputation = parallel;
    this->nbCore = nbCore;
    std::list<double> stploc;
    std::list<double> nrmloc;
    std::list<double> cstloc;
    std::list<double> tmeloc;
    std::list<int> hddnloc;
    double timeExc = get_wall_time();
    int displaySizeit = (int) std::log10(nMaxIter);
    displaySizeit ++;
    int displaySizehid = (int) std::log10(this->nPoints);
    displaySizehid ++;
    Eigen::VectorXd xStart = Eigen::VectorXd(this->nPoints);
    for(int i = 0; i < nPoints; i++){
        xStart[i] = 0.;
    }
                    
    heavyBall hBall = heavyBall(-1.0,gradTol,nMaxIter);
    uint i = 0;
    derivInfo current = this->oracle(xStart);
    derivInfo tmp;
    Eigen::VectorXd yk,xk,xkm1;
    xk = xStart;
    double cf, gNorm;
    gNorm = current.gradient.norm();
    cf = current.costFunction;
    int bestTau;

    nrmloc.push_back(gNorm);
    cstloc.push_back(cf);
    hddnloc.push_back(current.nbHidden);
    double timer;

    while(i< hBall.iMax() &&  gNorm > hBall.gradN()){
        timer = get_wall_time();
       if(i ==0){
            xkm1 = xk;
            bestTau = computeRanking(cf,gNorm,current.nbHidden,xk,current.gradient,step,nStep,scoring);
            xk = xk+step[bestTau]*current.gradient;
        }
        else{
            yk = xk + ((double) i - 1)/((double) i + 2)*(xk-xkm1);
            xkm1 = xk;
            tmp = this->oracle(yk);
            bestTau = computeRanking(current.costFunction,current.gradient.norm(),current.nbHidden,yk,tmp.gradient,step,nStep,scoring);
            xk = yk + step[bestTau]*tmp.gradient; 
        }
        current = this->oracle(xk);
        std::cout.precision(5);
        cf = current.costFunction;
        gNorm = current.gradient.norm();
        stploc.push_back(step[bestTau]);
        nrmloc.push_back(gNorm);
        cstloc.push_back(cf);
        hddnloc.push_back(current.nbHidden);
        if( i % 100 == 0){
            std::cout<<"i = "<<std::setfill('0')<<std::setw(displaySizeit)<<i<<" #Hidden = "<<current.ezHidden<<std::scientific<<" f = "<<cf<<" ||∇f|| = "<<gNorm<<" tau = "<<step[bestTau]<<std::endl;
        }
        i++;
        tmeloc.push_back(get_wall_time()-timer);
    }
    timeExc = get_wall_time() -timeExc;
    std::string cv = (i == hBall.iMax())? "✘" : "✓";
    std::cout<<"convergence "<<cv<<" Total time "<<std::fixed<<timeExc<<"sec Time per it "<<std::scientific<<timeExc/((double) i)<<"sec"<<std::endl;
    this->dataConvergence.stp = stploc;
    this->dataConvergence.nrm = nrmloc;
    this->dataConvergence.cst = cstloc;
    this->dataConvergence.tme = tmeloc;
    this->dataConvergence.hddn = hddnloc;
    this->dataConvergence.nbtt = ((int) i == nMaxIter)? -1 : (int) i;
    this->dataConvergence.wllTme = timeExc;
}

void OTproblem::perfomOptimPsiPolyakResearch(int lenAlpha,double* alpha,double gradTol,int nMaxIter,bool parallel,int nbCore){
    // introduction to optimisation Борис Поляк page 72 eqn(24)
    this->parallelComputation = parallel;
    this->nbCore = nbCore;
    std::list<double> stploc;
    std::list<double> nrmloc;
    std::list<double> cstloc;
    std::list<double> tmeloc;
    std::list<int> hddnloc;
    double timeExc = get_wall_time();
    int displaySizeit = (int) std::log10(nMaxIter);
    displaySizeit ++;
    int displaySizehid = (int) std::log10(this->nPoints);
    displaySizehid ++;
    Eigen::VectorXd xStart = Eigen::VectorXd(this->nPoints);
    for(int i = 0; i < nPoints; i++){
        xStart[i] = 0.;
    }
                    
    heavyBall hBall = heavyBall(-1.0,gradTol,nMaxIter);
    uint i = 0;
    derivInfo current = this->oracle(xStart);
    derivInfo tmp;
    Eigen::VectorXd xk,pk,pkm1,rk,rkm1;
    xk = xStart;
    double cf, gNorm,beta;
    gNorm = current.gradient.norm();
    cf = current.costFunction;
    nrmloc.push_back(gNorm);
    cstloc.push_back(cf);
    hddnloc.push_back(current.nbHidden);
    double timer;
    int alphaStar;
    while(i< hBall.iMax() &&  gNorm > hBall.gradN()){
        timer = get_wall_time();
       if(i ==0){
            beta = 0;
            rk = current.gradient;
            pk = rk;
            alphaStar = computeAlphaStar(cf,xk,pk,alpha,lenAlpha);
            xk = xk+alpha[alphaStar]*pk;
        }
        else{
            rk = current.gradient;
            beta = rk.dot(rk - rkm1)/rkm1.squaredNorm();
            pk = rk + beta*pkm1;
            alphaStar = computeAlphaStar(cf,xk,pk,alpha,lenAlpha);
            xk = xk+alpha[alphaStar]*pk;
        }
        pkm1 = pk; rkm1 = rk;
        current = this->oracle(xk);
        std::cout.precision(5);
        cf = current.costFunction;
        gNorm = current.gradient.norm();
        stploc.push_back(alpha[alphaStar]);
        nrmloc.push_back(gNorm);
        cstloc.push_back(cf);
        hddnloc.push_back(current.nbHidden);
        if( i % 100 ==0){
            std::cout<<"i = "<<std::setfill('0')<<std::setw(displaySizeit)<<i<<" #Hidden = "<<current.ezHidden<<std::scientific<<" f = "<<cf<<" ||∇f|| = "<<gNorm<<" tau = "<<alpha[alphaStar]<<std::endl;
        }
        i++;
        tmeloc.push_back(get_wall_time()-timer);
    }
    timeExc = get_wall_time() -timeExc;
    std::string cv = (i == hBall.iMax())? "✘" : "✓";
    std::cout<<"convergence "<<cv<<" Total time "<<std::fixed<<timeExc<<"sec Time per it "<<std::scientific<<timeExc/((double) i)<<"sec"<<std::endl;
    this->dataConvergence.stp = stploc;
    this->dataConvergence.nrm = nrmloc;
    this->dataConvergence.cst = cstloc;
    this->dataConvergence.tme = tmeloc;
    this->dataConvergence.hddn = hddnloc;
    this->dataConvergence.nbtt = ((int) i == nMaxIter)? -1 : (int) i;
    this->dataConvergence.wllTme = timeExc;
}

void OTproblem::perfomOptimPsiBFGSUnleashDaBeastResearch(double gradTol,int nMaxIter,double w1, double w2, int wMaxIter,int memSize,double leven,bool parallel,int nbCore){
    this->parallelComputation = parallel;
    this->nbCore = nbCore;
    std::list<double> stploc;
    std::list<double> nrmloc;
    std::list<double> cstloc;
    std::list<double> tmeloc;
    std::list<int> hddnloc;
    double timeExc = get_wall_time();
    int displaySizeit = (int) std::log10(nMaxIter);
    displaySizeit ++;
    int displaySizehid = (int) std::log10(this->nPoints);
    displaySizehid ++;
    Eigen::VectorXd xStart = Eigen::VectorXd(this->nPoints);
    for(int i = 0; i < nPoints; i++){
        xStart[i] = 0.;
    }
    int numberWithoutHidden=0;
    BFGS bf = BFGS((uint)memSize,(uint)this->nPoints);
    int i = 0;
    Eigen::VectorXd d;
    derivInfo current = this->oracle(xStart);
    derivInfo next;
    nrmloc.push_back(current.gradient.norm());
    cstloc.push_back(current.costFunction);
    hddnloc.push_back(current.nbHidden);
    Eigen::VectorXd x = xStart;
    double g;
    double gradGdotd;
    double s, sm, sp;
    int iterWolfe;
    double timer;
    while(i< nMaxIter && current.gradient.norm() > gradTol){
        timer = get_wall_time();
        // compute step s with wolfe conditions
        if ( numberWithoutHidden > 1){
            std::cout<<"Yo using Newton info"<<std::endl;
            matCOO coucou = this->Lag.Intersex(this->pol);
            Eigen::SparseMatrix<double> Hess =  coucou.toSparseEigen();
            Eigen::SparseMatrix<double> I(coucou.getDimRow(),coucou.getDimCol());
            I.setIdentity();
            I *= leven;
            try{
                Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> chol(Hess-I);
                d = chol.solve(-current.gradient);
            } catch (...){
                std::cout<<"some shit goes wrong"<<std::endl;
                d = bf.findDirection(current.gradient);    
            }
            

        }
        else{
            d = bf.findDirection(current.gradient);
        }
        
        s  = 1.0;
        sm =  0.0;
        sp = std::numeric_limits<double>::infinity();
        next = oracle(x+s*d);
        g  = current.costFunction;
        gradGdotd = current.gradient.dot(d);
        iterWolfe = 0;
        while( ((next.costFunction < g + w1*s*gradGdotd) || (next.gradient.dot(d) > w2*gradGdotd)) && iterWolfe < wMaxIter){
            if(next.costFunction < g + w1*s*gradGdotd){
                sp = s;
                s = (sm+sp)/2.;
            }
            else{
                sm=s;
                if(sp == std::numeric_limits<double>::infinity()){
                    s = 2.*s;
                }
                else{
                    s = (sm+sp)/2.;
                }
            }
            next = oracle(x + s*d);
            iterWolfe++;
        }
        std::string decreasecf   = (next.costFunction > g + w1*s*gradGdotd)? blue+"✓"+reset : red+"✘"+reset;
        std::string decreasegrad = (next.gradient.dot(d) < w2*gradGdotd)   ? blue+"✓"+reset : red+"✘"+reset;
        std::cout.precision(5);
        std::cout<<"i = "<<std::setfill('0')<<std::setw(displaySizeit)<<i<<" #Hidden = "<<current.ezHidden<<std::scientific<<" f = "<<current.costFunction<<" ||∇f|| = "<<current.gradient.norm();
        std::cout.precision(2);
        std::cout<<" s = "<<s<<" ↓cf : "<<decreasecf<<" ↓∇f "<<decreasegrad<<std::endl;
        bf.addDirection(s,d,next.gradient,current.gradient);
        
        current.costFunction = next.costFunction;
        current.gradient     = next.gradient;
        current.nbHidden     = next.nbHidden;
        current.ezHidden     = next.ezHidden;
        nrmloc.push_back(current.gradient.norm());
        cstloc.push_back(current.costFunction);
        hddnloc.push_back(current.nbHidden);
        stploc.push_back(s);
        x = x + s*d;
        i++;
        if (current.nbHidden == 0){
            numberWithoutHidden ++;
        }
        else{
            numberWithoutHidden = 0;
        }

        tmeloc.push_back(get_wall_time()-timer);
    }
    timeExc = get_wall_time() -timeExc;
    std::string cv = (i == nMaxIter)? "✘" : "✓";
    std::cout<<"convergence "<<cv<<" Total time "<<std::fixed<<timeExc<<"sec Time per it "<<std::scientific<<timeExc/((double) i)<<"sec"<<std::endl;
    this->optimizedPsi.phiOpt = x;
    this->optimizedPsi.costFunction = current.costFunction;
    this->optimizedPsi.wallTime = timeExc;
    this->optimizedPsi.nbIter = (i == nMaxIter)? -1 : i;
    this->dataConvergence.stp = stploc;
    this->dataConvergence.nrm = nrmloc;
    this->dataConvergence.cst = cstloc;
    this->dataConvergence.tme = tmeloc;
    this->dataConvergence.hddn = hddnloc;
    this->dataConvergence.nbtt = ((int) i == nMaxIter)? -1 : (int) i;
    this->dataConvergence.wllTme = timeExc;

}

void OTproblem::perfomOptimPsiHessLM(double gradTol,int nMaxIter,double w1, double w2, int wMaxIter,int memSize,double leven,bool parallel,int nbCore){
    this->parallelComputation = parallel;
    this->nbCore = nbCore;
    double c;
    std::list<double> stploc;
    std::list<double> nrmloc;
    std::list<double> cstloc;
    std::list<double> tmeloc;
    std::list<int> hddnloc;
    double timeExc = get_wall_time();
    int displaySizeit = (int) std::log10(nMaxIter);
    displaySizeit ++;
    int displaySizehid = (int) std::log10(this->nPoints);
    displaySizehid ++;
    Eigen::VectorXd xStart = Eigen::VectorXd(this->nPoints);
    for(int i = 0; i < nPoints; i++){
        xStart[i] = 0.;
    }
    int numberWithoutHidden=0;
    BFGS bf = BFGS((uint)memSize,(uint)this->nPoints);
    int i = 0;
    Eigen::VectorXd d,Dir;
    derivInfo current = this->oracle(xStart);
    derivInfo next;
    matCOO coucou = this->Lag.Intersex(this->pol);
    Eigen::SparseMatrix<double> Hess,Hess2;
    Hess =  coucou.toSparseEigen();
    Eigen::VectorXd diag = Hess.diagonal();
    c = 1e-2*diag.minCoeff();
    
    Eigen::SparseMatrix<double> I(coucou.getDimRow(),coucou.getDimCol());
    I.setIdentity();

    nrmloc.push_back(current.gradient.norm());
    cstloc.push_back(current.costFunction);
    hddnloc.push_back(current.nbHidden);
    Eigen::VectorXd x = xStart;
    double g;
    double predicted_increase;
    double s;
    double timer;
    int iterWolfe;
    while( (current.gradient.norm() > gradTol || current.nbHidden > 0) && i< nMaxIter ){
        timer = get_wall_time();
        c /= 1.7;
        Hess2 = Hess + c*I;
        try{
            Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> chol(Hess2);    
            Dir = chol.solve(-current.gradient);
        } catch (...){
            std::cout<<"some shit goes wrong"<<std::endl;
        }
        next = oracle(x+Dir);
        predicted_increase = current.gradient.dot(Dir);
        iterWolfe = 0;       
        while( ( ((next.gradient.norm()  > w2*current.gradient.norm()) && (next.costFunction < current.costFunction + w1*predicted_increase)) || next.nbHidden > current.nbHidden) && iterWolfe < wMaxIter ) {
            iterWolfe += 1;
            c *= 3.1;
            Hess2 = Hess + c*I;
            try{
                Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> chol(Hess2);    
                Dir = chol.solve(-current.gradient);
            } catch (...){
                std::cout<<"some shit goes wrong"<<std::endl;
            }   
            next = oracle(x+Dir);
            if (c < -1e15) {
                std::cout<<"godamn bouleshit Newton Direction"<<std::endl;
                Dir = 1e-3*current.gradient;
                next = oracle(x+Dir);
                c = -1e3;
                break;

            }            
        }

        g  = current.costFunction;
        
        s = Dir.norm();
        std::string decreasecf   = (next.costFunction > g + w1*predicted_increase)? blue+"✓"+reset : red+"✘"+reset;
        std::string decreasegrad = (next.gradient.dot(Dir) < w2*predicted_increase)   ? blue+"✓"+reset : red+"✘"+reset;
        std::cout.precision(5);
        std::cout<<"i = "<<std::setfill('0')<<std::setw(displaySizeit)<<i<<" #Hidden = "<<current.ezHidden<<std::scientific<<" f = "<<current.costFunction<<" ||∇f|| = "<<current.gradient.norm()<<" wolfe iter = "<<iterWolfe;
        std::cout.precision(2);
        std::cout<<" c = "<<c<<" ↓cf : "<<decreasecf<<" ↓∇f "<<decreasegrad<<std::endl;
        
        current.costFunction = next.costFunction;
        current.gradient     = next.gradient;
        current.nbHidden     = next.nbHidden;
        current.ezHidden     = next.ezHidden;
        nrmloc.push_back(current.gradient.norm());
        cstloc.push_back(current.costFunction);
        hddnloc.push_back(current.nbHidden);
        stploc.push_back(s);
        x = x + Dir;
        i++;
        if (current.nbHidden == 0){
            numberWithoutHidden ++;
        }
        else{
            numberWithoutHidden = 0;
        }

        tmeloc.push_back(get_wall_time()-timer);
    }
    timeExc = get_wall_time() -timeExc;
    std::string cv = (i == nMaxIter)? "✘" : "✓";
    std::cout<<"convergence "<<cv<<" Total time "<<std::fixed<<timeExc<<"sec Time per it "<<std::scientific<<timeExc/((double) i)<<"sec"<<std::endl;
    this->optimizedPsi.phiOpt = x;
    this->optimizedPsi.costFunction = current.costFunction;
    this->optimizedPsi.wallTime = timeExc;
    this->optimizedPsi.nbIter = (i == nMaxIter)? -1 : i;
    this->dataConvergence.stp = stploc;
    this->dataConvergence.nrm = nrmloc;
    this->dataConvergence.cst = cstloc;
    this->dataConvergence.tme = tmeloc;
    this->dataConvergence.hddn = hddnloc;
    this->dataConvergence.nbtt = ((int) i == nMaxIter)? -1 : (int) i;
    this->dataConvergence.wllTme = timeExc;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////  ____                   _           __                  _   _                 ///////////////
/////////////// / ___|_ __ ___  ___  __| |_   _    / _|_   _ _ __   ___| |_(_) ___  _ __  ___ ///////////////
///////////////| |  _| '__/ _ \/ _ \/ _` | | | |  | |_| | | | '_ \ / __| __| |/ _ \| '_ \/ __|///////////////
///////////////| |_| | | |  __/  __/ (_| | |_| |  |  _| |_| | | | | (__| |_| | (_) | | | \__ \///////////////
 ///////////////\____|_|  \___|\___|\__,_|\__, |  |_|  \__,_|_| |_|\___|\__|_|\___/|_| |_|___////////////////
////////////////////////////////////////  |___/  ////////////////////////////////////////////////////////////


int OTproblem::computeRanking(double cf,double gNorm, int nbHid,const Eigen::VectorXd& xk, const Eigen::VectorXd& direction, double * tau,int lenTau, double * score){
    // il faut paralleliser cette section avec des threads
    int indexBest = 0;
    double BestScore = -100;
    double localScore,cfLocal,gNormLocal;
    derivInfo tmp;
    std::deque<int>    nbHiddenList;
    std::deque<double> cfList;
    std::deque<double> gnList;
    std::deque<double> scoreList;
    for(int i = 0; i < lenTau; i ++){
        localScore = 0;
        tmp = this->oracle(xk+tau[i]*direction);
        cfLocal = tmp.costFunction;
        gNormLocal = tmp.gradient.norm();
        nbHiddenList.push_back(tmp.nbHidden); cfList.push_back(cfLocal); gnList.push_back(gNormLocal);
        localScore += cfLocal > cf ? score[0] : 0.;
        localScore += gNormLocal < gNorm ? score[1] : 0.;
        if (tmp.nbHidden < nbHid){
            localScore += score[2] < 0 ? -score[2] : score[2]* ((double) nbHid - tmp.nbHidden);
        }
        else{
            if(score[2] > 0){
                localScore += score[2]* ((double) nbHid - tmp.nbHidden);
            } 
        }
        scoreList.push_back(localScore);
        if ( localScore > BestScore ){
            indexBest = i;
            BestScore = localScore;

        }
    }
    
    std::list<int> sameScore;
    sameScore.push_back(indexBest);
    int nb = 0;
    for( const auto & s : scoreList){
        if( std::abs(s - BestScore) < 1e-12){
            sameScore.push_back(nb);
        }
        nb++;
    }
    if( sameScore.size() > 1 ){
        // there is multiple tau with the same value choosing the one smaller GN
        // we could potentially use a non-binary relation for ranking ascending cf
        // and decreasing gradient norm
        for( auto it  = ++sameScore.begin(); it != sameScore.end(); ++it){
            if( gnList[indexBest] > gnList[*it]){
                indexBest = *it;
            }
        }
    }

    // std::cout.precision(5);
    // for( int inde = 0; inde < lenTau;inde++){
    //     std::cout<<"i = "<<inde<<std::scientific<<" tau = "<<tau[inde]<<" cf ↑ "<<cfList[inde]-cf<<" gradient ↓ "<<gNorm-gnList[inde]<<" Hidden ⇵ "<<nbHiddenList[inde]<<"/"<<nbHid<<" score "<<scoreList[inde]<<std::endl;
    // }
    // std::cout<<"Best is : "<<indexBest<<std::endl;

    return indexBest;
}



int OTproblem::computeRankingCont(double cf,double gNorm, int nbHid,const Eigen::VectorXd& xk, const Eigen::VectorXd& direction, double * tau,int lenTau, double * score){
    int indexBest = 0;
    double BestScore = -100;
    double localScore,cfLocal,gNormLocal;
    derivInfo tmp;
    std::deque<int>    nbHiddenList;
    std::deque<double> cfList;
    std::deque<double> gnList;
    std::deque<double> scoreList;
    for(int i = 0; i < lenTau; i ++){
        localScore = 0;
        tmp = this->oracle(xk+tau[i]*direction);
        cfLocal = tmp.costFunction;
        gNormLocal = tmp.gradient.norm();
        nbHiddenList.push_back(tmp.nbHidden); cfList.push_back(cfLocal); gnList.push_back(gNormLocal);
        localScore += (cfLocal - cf)/std::abs(cfLocal)*score[0];
        localScore += (gNorm-gNormLocal)/std::abs(gNormLocal)*score[1];        
        if (tmp.nbHidden < nbHid){
            localScore += score[2] < 0 ? -score[2] : score[2]* ((double) nbHid - tmp.nbHidden);
        }
        else{
            if(score[2] > 0){
                localScore += score[2]* ((double) nbHid - tmp.nbHidden);
            } 
        }
        scoreList.push_back(localScore);
        if ( localScore > BestScore ){
            indexBest = i;
            BestScore = localScore;

        }
    }
    
    std::list<int> sameScore;
    sameScore.push_back(indexBest);
    int nb = 0;
    for( const auto & s : scoreList){
        if( std::abs(s - BestScore) < 1e-12){
            sameScore.push_back(nb);
        }
        nb++;
    }
    if( sameScore.size() > 1 ){
        // there is multiple tau with the same value choosing the one smaller GN
        // we could potentially use a non-binary relation for ranking ascending cf
        // and decreasing gradient norm
        for( auto it  = ++sameScore.begin(); it != sameScore.end(); ++it){
            if( gnList[indexBest] > gnList[*it]){
                indexBest = *it;
            }
        }
    }

    // std::cout.precision(5);
    // for( int inde = 0; inde < lenTau;inde++){
    //     std::cout<<"i = "<<inde<<std::scientific<<" tau = "<<tau[inde]<<" cf ↑ "<<cfList[inde]-cf<<" gradient ↓ "<<gNorm-gnList[inde]<<" Hidden ⇵ "<<nbHiddenList[inde]<<"/"<<nbHid<<" score "<<scoreList[inde]<<std::endl;
    // }
    // std::cout<<"Best is : "<<indexBest<<std::endl;

    return indexBest;
}

int OTproblem::computeAlphaStar(double cf ,const Eigen::VectorXd& xk, const Eigen::VectorXd& direction, double * alpha,int lenAplha){
    // il faut paralleliser cette section avec des threads
    derivInfo tmp = this->oracle(xk+alpha[0]*direction);
    int indexBest = 0;
    double bestCostFunction = tmp.costFunction;
    for( int i = 1; i < lenAplha; i++){
        tmp = this->oracle(xk+alpha[i]*direction);
        if( bestCostFunction < tmp.costFunction){
            bestCostFunction = tmp.costFunction;
            indexBest = i;
        }
    }
    if( bestCostFunction < cf){
        std::cout<<this->red<<"Watch out no improvement of the cost function found in this direction"<<this->reset<<std::endl;
    }
    return indexBest;
}




//////////////////////////////////////////////////////////////////
/////////  ____ _____ _____ _____ _____ ____  ____  //////////////
///////// / ___| ____|_   _|_   _| ____|  _ \/ ___| //////////////
/////////| |  _|  _|   | |   | | |  _| | |_) \___ \ //////////////
/////////| |_| | |___  | |   | | | |___|  _ < ___) |//////////////
///////// \____|_____| |_|   |_| |_____|_| \_\____/ //////////////
//////////////////////////////////////////////////////////////////                                           



void OTproblem::getOptimizedPsi(int nMass,double* mass){
    if(nMass != this->nPoints){
        throw std::invalid_argument( "Dimension mismatch In getOptimizedPsi");
    }
    for(int i = 0; i < nMass; i++){
        mass[i] = this->optimizedPsi.phiOpt[i];
    }
}
double OTproblem::getOptimizationtime(){
    return this->optimizedPsi.wallTime;
}
double OTproblem::getOptimizationCF(){
    return this->optimizedPsi.costFunction;
}
int OTproblem::getOptimizationCV(){
    return this->optimizedPsi.nbIter;
}
std::list<double> OTproblem::getContributorLine(int i){
    int countPuff = 0;    
    for( auto & l : this->pol.lines ){
        if( i == countPuff )
            return l.infoBarycenter;
        countPuff ++;
    }
    return std::list<double>(0);   
}
