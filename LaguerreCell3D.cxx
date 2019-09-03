#include <stdexcept>
#include "LaguerreCell3D.h"
#include <assert.h>

#include <sys/time.h> 

#include <thread>

// MISC f
double square(Point p1){return std::pow(p1.x(),2.0) + std::pow(p1.y(),2.0) + std::pow(p1.z(),2.0);};
double dot(Point p1, Point p2){return p1.x()*p2.x()+ p1.y()*p2.y() + p1.z()*p2.z();};

// Flush functions
void LaguerreCell3D::flush(){    
    for (int i=0;i<this->nPoints; ++i){
        primalPtList[i].flush();      
    }
    this->OriginInstanciated=false;
};

// code
LaguerreCell3D::LaguerreCell3D(){
    this->PsiInstanciated=false;
    this->CoordInstanciated=false;
    this->OriginInstanciated=false;
}

LaguerreCell3D::LaguerreCell3D(const LaguerreCell3D& other){
    this->nPoints = other.nPoints;
    this->origin  = other.origin;
    this->originLocation = other.originLocation; 
    this->primalPtList = other.primalPtList;
    this->CoordInstanciated = true;
    this->PsiInstanciated = true;
    this->OriginInstanciated = true;
}

void LaguerreCell3D::setCoordinates(double* positionX,double* positionY,double* positionZ,int nX){
    this->primalPtList.clear();
	this->nPoints = nX;
    for (int i = 0; i < this->nPoints; ++i){
        this->primalPtList.push_back(PrimalPoint(Point(positionX[i],positionY[i],positionZ[i]),i));
    }
    this->CoordInstanciated=true;
    this->flush();    
    this->PsiInstanciated=false;
    this->OriginInstanciated=false;
}

void LaguerreCell3D::ComputeTesselation3D(double* weight,int nW){
    this->flush();
    if (!this->CoordInstanciated){
		throw std::invalid_argument( "No coordinates instanciated");
    }
    if(this->nPoints != nW){
        throw std::invalid_argument( "Number of points and psi does not match" );
    }
	for (int i = 0; i < this->nPoints; ++i){
        this->primalPtList[i].w(weight[i]);
    }
    
    std::vector<std::pair<Weighted_point,int> > points_tmp;
    for (const auto& PrimalPt :  this->primalPtList){
        points_tmp.push_back(std::make_pair(Weighted_point(PrimalPt.point,PrimalPt.w()),PrimalPt.id));
    }
    T.clear();
    T.insert(points_tmp.begin(),points_tmp.end());
    T.infinite_vertex()->info() = -1;
    // timeExc = get_wall_time() -timeExc; 
    // std::cout<<"from CGAL : "<<timeExc<<std::endl;
    // timeExc = get_wall_time();
    if (T.dimension()==2){
        this->GetdualInfo2d();
    }
    else {
        this->GetdualInfo3d();
    }
    // timeExc = get_wall_time() -timeExc; 
    // std::cout<<"from GetDualInfo : "<<timeExc<<std::endl;
    this->PsiInstanciated=true;
    this->OriginInstanciated=false;
    this->SetOrigin();
}
void LaguerreCell3D::GetdualInfo3d(){
    for (Cell_iterator cit = T.cells_begin(); cit != T.cells_end(); cit++){
        std::list<int> finiteVertices(0);
        for (int j = 0; j < 4; ++j){
            if (!T.is_infinite(cit->vertex(j))){
                finiteVertices.push_back(cit->vertex(j)->info());
            }
        }
        for (auto j : finiteVertices){
            for( auto testj : finiteVertices){
                if( testj != j) this->primalPtList[j].adjPrimalPoints.insert(testj);
            }           
        }
        
    }
}
void LaguerreCell3D::GetdualInfo2d(){
    Point PtTmp;
    for (Finite_facets_iterator fit = T.finite_facets_begin(); fit != T.finite_facets_end(); fit++){
        std::list<int> finiteVertices;
        for (int j = 0; j < 3; ++j){
            finiteVertices.push_back(fit->first->vertex( (fit->second+j+1)%4)->info());
        }
        for (auto j : finiteVertices){
            for( auto testj : finiteVertices){
                if( testj != j) this->primalPtList[j].adjPrimalPoints.insert(testj);
                
            }           
        }
    }
}

double LaguerreCell3D::ComputeEnergy(const Point& P,const PrimalPoint& Pi){
    // Compute || P-P_i||^2-w_i if P=(X,Y,Z) and (P_i,w_i) is the weighted point indice i
    if (!this->CoordInstanciated || !this->PsiInstanciated){
        throw std::invalid_argument( "Coordinates and Psi has to be instanciated before using this function");
    }
    return std::pow(P.x()-Pi.x(),2.0)+ std::pow(P.y()-Pi.y(),2.0)+std::pow(P.z()-Pi.z(),2.0)-Pi.w();
    // return (P-Pi.point).squared()-Pi.w()
}

bool LaguerreCell3D::IsPointInCell(const Point& P, const PrimalPoint& Pi){
    double energ=ComputeEnergy(P,Pi);
    bool isInside=true;
    for (auto& j  : Pi.adjPrimalPoints){
        double energ2=ComputeEnergy(P,this->primalPtList[j]);
        if (energ2 <energ){
            isInside=false;
        }
    }
    return isInside;
}


void LaguerreCell3D::SetOrigin(){
    int ind=0;
    double max=this->primalPtList[0].w();
    for (auto & P: this->primalPtList){
        if (P.w()>max){
            max=P.w();
            ind=P.id;
        }
    }
    this->origin=this->primalPtList[ind].point;
    this->originLocation=ind;
    this->OriginInstanciated=true;
    if (!IsPointInCell(this->origin,this->primalPtList[this->originLocation])){
        std::cout<<"Could NOT FIND ORIGIN "<<ind<<" "<<max<<std::endl;
    }
}

int LaguerreCell3D::Locate(const Point& P,const Point& surePoint,int sureLocation){
    double actualTime=0.;
    int actualLoc=sureLocation;
    GenericResult r;
    int previousLoc = -1;
    while (!(actualTime==1.)){
        r=FindNextIntersection(surePoint,P, actualTime,actualLoc,previousLoc); 
        previousLoc = actualLoc;    
        actualLoc=r.location;
        actualTime=r.time;
    }
    return r.location;
}

GenericResult LaguerreCell3D::FindNextIntersection(const Point& P0,const Point& P1, double actualTime, int actualLoc,int previousLoc){
    GenericResult r;
    bool hasExitTime=false;
    r.time=1.;
    r.location=actualLoc;
    double exitTime;
    PrimalPoint& Pi=primalPtList[actualLoc];
    for (auto& j : this->primalPtList[actualLoc].adjPrimalPoints){
        PrimalPoint& Pj=primalPtList[j];
        if( j == previousLoc) continue;
        double fact=2*( (P0.x()-P1.x())*(Pi.x()-Pj.x())+
                        (P0.y()-P1.y())*(Pi.y()-Pj.y())+
                        (P0.z()-P1.z())*(Pi.z()-Pj.z()));
                        // REMPLACER PAR fact=2*dot(P0-P1,Pi.point-Pj.point)
        if (std::abs(fact) >1.e-10){
            double tmp=Pi.w()-Pj.w();
//            tmp-=std::pow(Pi.x(),2.0)+std::pow(Pi.y(),2.0)+std::pow(Pi.z(),2.0);
            tmp-=Pi.x()*Pi.x()+Pi.y()*Pi.y()+Pi.z()*Pi.z();
//            tmp+=std::pow(Pj.x(),2.0)+std::pow(Pj.y(),2.0)+std::pow(Pj.z(),2.0); 
            tmp+=Pj.x()*Pj.x()+Pj.y()*Pj.y()+Pj.z()*Pj.z(); 
            tmp+=2*( P0.x()*(Pi.x()-Pj.x())+
                     P0.y()*(Pi.y()-Pj.y())+
                     P0.z()*(Pi.z()-Pj.z()));
            // REMPLACER PAR TMP= ComputeEnergy(P0,Pj)-ComputeEnergy(P0,Pi)
            exitTime=tmp/fact;
            if (exitTime > actualTime){
                hasExitTime=true;
                if (exitTime <r.time) {
                    r.time=exitTime;
                    r.location=j;
                }
            }
        }
    }
    if (!hasExitTime){
        r.time=1.;
    }
    if (r.time==1.){
        if (!IsPointInCell(P1, Pi)){
            std::cout<<"ENDING POINT IS NOT IN THE POLYGON AND SHOULD BE. hasexited ?"<<hasExitTime<<std::endl;
            throw std::invalid_argument("dd");
        }
    }
    return r;


}

void LaguerreCell3D::Intersect(Polyline& polyline){
    GenericResult r;    
    polyline.flush();
    int actualLoc=this->Locate(polyline.lines.front().startPoint,this->origin,this->originLocation);
    Point * SurePoint =&polyline.lines.front().startPoint;
    Point Pbarycenter;
    double lengthLine;
    double localMass,localCost;
    double gaussTimeLeft, gaussTimeRight;
    Point gaussPointLeft, gaussPointRight;
    PrimalPoint * currentPrimalPoint;
    int previousLoc;
    for (auto l=polyline.lines.begin();l!=polyline.lines.end();l++){
        actualLoc=this->Locate(l->startPoint,*(SurePoint),actualLoc);
        previousLoc = -1;   
        double actualTime=0.;
        std::list <int> intersectionLocation;
        std::list <double> intersectionTime;
        intersectionLocation.push_back(actualLoc);
        intersectionTime.push_back(actualTime);     
        lengthLine = l->lengthLine;

        
        gaussTimeLeft   = (0.5 - 0.5/std::sqrt(3));
        gaussTimeRight  = (0.5 + 0.5/std::sqrt(3));
        gaussPointLeft  = l->getPointFromTime(gaussTimeLeft);
        gaussPointRight = l->getPointFromTime(gaussTimeRight);
        
        while (!(actualTime==1.)){
            r=FindNextIntersection(l->startPoint,l->endPoint, actualTime,actualLoc,previousLoc);
            currentPrimalPoint = &primalPtList[actualLoc];
            localMass = (r.time-actualTime)*l->weight*lengthLine;
            Pbarycenter = l->getPointFromTime((r.time+actualTime)*.5);
            currentPrimalPoint->barycenter[0] += (Pbarycenter.x())*localMass;
            currentPrimalPoint->barycenter[1] += (Pbarycenter.y())*localMass;
            currentPrimalPoint->barycenter[2] += (Pbarycenter.z())*localMass;
              
            gaussTimeLeft = actualTime + (0.5 - 0.5/std::sqrt(3))*(r.time-actualTime);
            gaussTimeRight = actualTime + (0.5 + 0.5/std::sqrt(3))*(r.time-actualTime);
            gaussPointLeft  = l->getPointFromTime(gaussTimeLeft);
            gaussPointRight = l->getPointFromTime(gaussTimeRight);
            localCost = 0.5*(l->computeSquaredDistance(gaussTimeLeft,*currentPrimalPoint) 
                           + l->computeSquaredDistance(gaussTimeRight,*currentPrimalPoint))*localMass;
            
            currentPrimalPoint->mass += localMass;
            l->massSeen              -= localMass*currentPrimalPoint->w();
            currentPrimalPoint->cost += localCost;
            l->cost                  += localCost;
            
            previousLoc = actualLoc;
            actualLoc=r.location;
            actualTime=r.time;
            intersectionLocation.push_back(actualLoc);
            intersectionTime.push_back(actualTime);
        }

        int lastLoc = intersectionLocation.back();
        intersectionLocation.pop_back();        
        assert(lastLoc==intersectionLocation.back());
        SurePoint = &(l->endPoint);
        l->intersectionLocation.reserve(intersectionLocation.size());
        std::copy(std::begin(intersectionLocation), std::end(intersectionLocation), std::back_inserter(l->intersectionLocation));
        l->intersectionTime.reserve(intersectionTime.size());
        std::copy(std::begin(intersectionTime), std::end(intersectionTime), std::back_inserter(l->intersectionTime));

    }
    for( auto &primalPt : this->primalPtList){
        for( int i = 0; i < 3;i ++){
            primalPt.barycenter[i] /=  (primalPt.mass > 1e-10) ? primalPt.mass : 1 ;
        }        
    }
}

void LaguerreCell3D::Intersect2(Polyline& polyline){
    GenericResult r;    
    polyline.flush();
    int actualLoc=this->Locate(polyline.lines.front().startPoint,this->origin,this->originLocation);
    Point * SurePoint =&polyline.lines.front().startPoint;
    Point Pbarycenter;
    double lengthLine;
    double localMass,localCost;
    double gaussTimeLeft, gaussTimeRight;
    Point gaussPointLeft, gaussPointRight;
    PrimalPoint * currentPrimalPoint;
    int previousLoc;
    for (auto l=polyline.lines.begin();l!=polyline.lines.end();l++){
        actualLoc=this->Locate(l->startPoint,*(SurePoint),actualLoc);
        previousLoc = -1;   
        double actualTime=0.;
        std::list <int> intersectionLocation;
        std::list <double> intersectionTime;
        intersectionLocation.push_back(actualLoc);
        intersectionTime.push_back(actualTime);     
        lengthLine = l->lengthLine;

        
        gaussTimeLeft   = (0.5 - 0.5/std::sqrt(3));
        gaussTimeRight  = (0.5 + 0.5/std::sqrt(3));
        gaussPointLeft  = l->getPointFromTime(gaussTimeLeft);
        gaussPointRight = l->getPointFromTime(gaussTimeRight);
        
        while (!(actualTime==1.)){
            r=FindNextIntersection(l->startPoint,l->endPoint, actualTime,actualLoc,previousLoc);
            currentPrimalPoint = &primalPtList[actualLoc];
            localMass = (r.time-actualTime)*l->weight*lengthLine;
            Pbarycenter = l->getPointFromTime((r.time+actualTime)*.5);
            l->infoBarycenter.push_back((double) actualLoc);
            l->infoBarycenter.push_back(localMass);
            currentPrimalPoint->barycenter[0] += (Pbarycenter.x())*localMass;
            currentPrimalPoint->barycenter[1] += (Pbarycenter.y())*localMass;
            currentPrimalPoint->barycenter[2] += (Pbarycenter.z())*localMass;
              
            gaussTimeLeft = actualTime + (0.5 - 0.5/std::sqrt(3))*(r.time-actualTime);
            gaussTimeRight = actualTime + (0.5 + 0.5/std::sqrt(3))*(r.time-actualTime);
            gaussPointLeft  = l->getPointFromTime(gaussTimeLeft);
            gaussPointRight = l->getPointFromTime(gaussTimeRight);
            localCost = 0.5*(l->computeSquaredDistance(gaussTimeLeft,*currentPrimalPoint) 
                           + l->computeSquaredDistance(gaussTimeRight,*currentPrimalPoint))*localMass;
            
            currentPrimalPoint->mass += localMass;
            l->massSeen              -= localMass*currentPrimalPoint->w();
            currentPrimalPoint->cost += localCost;
            l->cost                  += localCost;
            
            previousLoc = actualLoc;
            actualLoc=r.location;
            actualTime=r.time;
            intersectionLocation.push_back(actualLoc);
            intersectionTime.push_back(actualTime);
        }

        int lastLoc = intersectionLocation.back();
        intersectionLocation.pop_back();        
        assert(lastLoc==intersectionLocation.back());
        SurePoint = &(l->endPoint);
        l->intersectionLocation.reserve(intersectionLocation.size());
        std::copy(std::begin(intersectionLocation), std::end(intersectionLocation), std::back_inserter(l->intersectionLocation));
        l->intersectionTime.reserve(intersectionTime.size());
        std::copy(std::begin(intersectionTime), std::end(intersectionTime), std::back_inserter(l->intersectionTime));

    }
    for( auto &primalPt : this->primalPtList){
        for( int i = 0; i < 3;i ++){
            primalPt.barycenter[i] /=  (primalPt.mass > 1e-10) ? primalPt.mass : 1 ;
        }        
    }
}


matCOO LaguerreCell3D::Intersex(Polyline& polyline){
    matCOO ret = matCOO(this->primalPtList.size());
    std::vector<double> dL(3);
    std::vector<double> dBorder(3);
    for( auto & l : polyline.lines){
        Point & Pe = l.endPoint;
        Point & Ps = l.startPoint;
        double lengthLine = l.lengthLine;
        dL[0]=(Pe.x()-Ps.x())/lengthLine;
        dL[1]=(Pe.y()-Ps.y())/lengthLine;
        dL[2]=(Pe.z()-Ps.z())/lengthLine;
        
        for(unsigned i = 0; i < l.intersectionLocation.size()-1; i++){
            PrimalPoint & current = this->primalPtList[l.intersectionLocation[i]];
            PrimalPoint & next = this->primalPtList[l.intersectionLocation[i+1]];
            dBorder[0] = next.x() - current.x();
            dBorder[1] = next.y() - current.y();
            dBorder[2] = next.z() - current.z();
            double dot = 0;
            for( int k = 0; k < 3; k++){
                dot += dBorder[k]*dL[k];
            }
            double value = .5*std::abs(l.weight/dot);
            ret.push(l.intersectionLocation[i]  ,l.intersectionLocation[i+1], value);
            ret.push(l.intersectionLocation[i+1],l.intersectionLocation[i]  , value);
            ret.push(l.intersectionLocation[i]  ,l.intersectionLocation[i]  ,-value);
            ret.push(l.intersectionLocation[i+1],l.intersectionLocation[i+1],-value);

        }
    }
    return ret;
}


void intersectionJob(LaguerreCell3D * L3,Polyline * polyline){
    GenericResult r;
    int actualLoc=L3->Locate(polyline->lines.front().startPoint,L3->origin,L3->originLocation);
    Point * SurePoint =&polyline->lines.front().startPoint;
    Point Pbarycenter;
    double lengthLine;
    double localMass,localCost;
    double gaussTimeLeft, gaussTimeRight;
    Point gaussPointLeft, gaussPointRight;
    PrimalPoint * currentPrimalPoint;
    int previousLoc;
    for (auto l=polyline->lines.begin();l!=polyline->lines.end();l++){
        actualLoc=L3->Locate(l->startPoint,*(SurePoint),actualLoc);
        previousLoc = -1;   
        double actualTime=0.;
        std::list <int> intersectionLocation;
        std::list <double> intersectionTime;
        intersectionLocation.push_back(actualLoc);
        intersectionTime.push_back(actualTime);     
        lengthLine = l->lengthLine;

        
        gaussTimeLeft   = (0.5 - 0.5/std::sqrt(3));
        gaussTimeRight  = (0.5 + 0.5/std::sqrt(3));
        gaussPointLeft  = l->getPointFromTime(gaussTimeLeft);
        gaussPointRight = l->getPointFromTime(gaussTimeRight);
        
        while (!(actualTime==1.)){
            r=L3->FindNextIntersection(l->startPoint,l->endPoint, actualTime,actualLoc,previousLoc);
            currentPrimalPoint = &L3->primalPtList[actualLoc];
            localMass = (r.time-actualTime)*l->weight*lengthLine;
              
            gaussTimeLeft = actualTime + (0.5 - 0.5/std::sqrt(3))*(r.time-actualTime);
            gaussTimeRight = actualTime + (0.5 + 0.5/std::sqrt(3))*(r.time-actualTime);
            gaussPointLeft  = l->getPointFromTime(gaussTimeLeft);
            gaussPointRight = l->getPointFromTime(gaussTimeRight);
            localCost = 0.5*(l->computeSquaredDistance(gaussTimeLeft,*currentPrimalPoint) 
                           + l->computeSquaredDistance(gaussTimeRight,*currentPrimalPoint))*localMass;
            
            currentPrimalPoint->mass += localMass;
            l->massSeen              -= localMass*currentPrimalPoint->w();
            currentPrimalPoint->cost += localCost;
            l->cost                  += localCost;
            
            previousLoc = actualLoc;
            actualLoc=r.location;
            actualTime=r.time;
            intersectionLocation.push_back(actualLoc);
            intersectionTime.push_back(actualTime);
        }

        int lastLoc = intersectionLocation.back();
        intersectionLocation.pop_back();        
        assert(lastLoc==intersectionLocation.back());
        SurePoint = &(l->endPoint);
        l->intersectionLocation.reserve(intersectionLocation.size());
        std::copy(std::begin(intersectionLocation), std::end(intersectionLocation), std::back_inserter(l->intersectionLocation));
        l->intersectionTime.reserve(intersectionTime.size());
        std::copy(std::begin(intersectionTime), std::end(intersectionTime), std::back_inserter(l->intersectionTime));

    }
}

void LaguerreCell3D::IntersectParallel(Polyline& polyline,int nbThread){
    polyline.flush();
    std::vector<LaguerreCell3D> classList(nbThread);
    std::vector<Polyline> polyList(nbThread);
    std::vector<std::thread> vectorThreads(nbThread);
    int n = polyline.lines.size();
    std::list<Line>::iterator deb;
    std::list<Line>::iterator fin;
    deb = polyline.lines.begin();
    fin = polyline.lines.begin();

    for(int m = 0; m < (n/nbThread); m++,fin++){}
    for(int i = 0; i < nbThread; i++){
        classList[i] = LaguerreCell3D((*this));
        // malheureusement pour les listes l'operation std::list<L>::iterator + int 
        // n'est pas acceptable, il faut donc construire ces listes a la main.
        // c'est seulement de bi directionnal iterator
        if( i != nbThread-1){
            // std::cout<<"Launching tread with istart "<<ideb<<", ifin "<<ifin<<std::endl;
            polyList[i] = Polyline(deb,fin);
            deb = fin;
            for(int m = 0; m < (n/nbThread); m++,fin++){}
        }
        else{
            // std::cout<<"Launching tread with istart "<<ideb<<", ifin "<<n-1<<std::endl;
            polyList[i] = Polyline(deb,polyline.lines.end());
        }
    }
    
    for(int t = 0; t < nbThread ; t++){
        vectorThreads[t] = std::thread(intersectionJob,&classList[t],&polyList[t]);
    }

    for(int t = 0; t < nbThread; t ++){
   		vectorThreads[t].join();
   	}
    
    for(int t = 0; t < nbThread; t ++){
        for(int i= 0; i < this->nPoints;i++){
   		    primalPtList[i].mass += classList[t].primalPtList[i].mass;
            primalPtList[i].cost += classList[t].primalPtList[i].cost;

        }
   	}
}


double Line::computeSquaredDistance(double t, const PrimalPoint & P){
    return std::pow((1-t)*this->startPoint.x() + t*this->endPoint.x() - P.x(),2)
          +std::pow((1-t)*this->startPoint.y() + t*this->endPoint.y() - P.y(),2)
          +std::pow((1-t)*this->startPoint.z() + t*this->endPoint.z() - P.z(),2);
}

void LaguerreCell3D::ComputePolyLineDer(Polyline& polyline){
    double t;
    std::vector<double> derivTStart, derivTend;
    derivTStart = std::vector<double>(6);
    derivTend   = std::vector<double>(6);
    double denominator, intePs, intePe, PmiddleDotxi;
    std::vector<double> tprime = std::vector<double>(6);
    std::vector<double> dL = std::vector<double>(6);
    
    for( auto & l : polyline.lines){
        Point & Pe = l.endPoint;
        Point & Ps = l.startPoint;
        double lengthLine = l.lengthLine;
        dL[0]=(Pe.x()-Ps.x())/lengthLine;
        dL[1]=(Pe.y()-Ps.y())/lengthLine;
        dL[2]=(Pe.z()-Ps.z())/lengthLine;
        double normSquaredPs = std::pow(Ps.x(),2) + std::pow(Ps.y(),2) + std::pow(Ps.z(),2);
        double normSquaredPe = std::pow(Pe.x(),2) + std::pow(Pe.y(),2) + std::pow(Pe.z(),2);
        double PedotPs = Ps.x()*Pe.x() + Ps.y()*Pe.y() + Ps.z()*Pe.z();
        // i
        for (int k=0;k<3;k++){
            l.costDerivStart[k] -= 1./3.*dL[k]*(normSquaredPe + normSquaredPs + PedotPs); //OK
            l.costDerivEnd[k]   += 1./3.*dL[k]*(normSquaredPe + normSquaredPs + PedotPs); //OK

        }
        
        // ii
        l.costDerivStart[0] += 1./3.*lengthLine*(2.*Ps.x() + Pe.x()); //OK
        l.costDerivStart[1] += 1./3.*lengthLine*(2.*Ps.y() + Pe.y()); //OK
        l.costDerivStart[2] += 1./3.*lengthLine*(2.*Ps.z() + Pe.z()); //OK

        l.costDerivEnd[0]   += 1./3.*lengthLine*(Ps.x() + 2*Pe.x()); //OK
        l.costDerivEnd[1]   += 1./3.*lengthLine*(Ps.y() + 2*Pe.y()); //OK
        l.costDerivEnd[2]   += 1./3.*lengthLine*(Ps.z() + 2*Pe.z()); //OK

        
        PrimalPoint & PPrevious =this->primalPtList[l.intersectionLocation[l.intersectionLocation.size()-1]];
        double squarePP = std::pow(PPrevious.x(),2.) + std::pow(PPrevious.y(),2.) + std::pow(PPrevious.z(),2.);
        for (int k=0;k<3;k++){
            l.massDerivStart[k]+=PPrevious.w()*dL[k];
            l.massDerivEnd[k]  -=PPrevious.w()*dL[k];
            l.costDerivStart[k] -=  squarePP*dL[k]; 
            l.costDerivEnd[k]   +=  squarePP*dL[k]; 
        }
        double ts,te;
        double diffTime;
        for(unsigned i = 0; i < l.intersectionLocation.size(); i++){
            ts = l.intersectionTime[i];
            te = l.intersectionTime[i+1];
            PrimalPoint & xi = this->primalPtList[l.intersectionLocation[i]];
            diffTime = te - ts;
            intePs = diffTime*(1.0-0.5*(te+ts));
            intePe = diffTime*(0.5*(te+ts));
        
            Point Pmiddle = l.getPointFromTime(0.5*(te+ts));
            PmiddleDotxi = Pmiddle.x()*xi.x() + Pmiddle.y()*xi.y() + Pmiddle.z()*xi.z();
            for (int k=0; k<3;k++){
                /*
                -2 dot(\| p_e - p_s\|) int_ts,te xt.xi
                */
                l.costDerivStart[k] += 2.*diffTime*dL[k]*PmiddleDotxi;
                l.costDerivEnd[k]   -= 2.*diffTime*dL[k]*PmiddleDotxi;
          
            }
            /*
                -2 \| p_e - p_s\| int_ts,te dot(xt).xi
            */
            l.costDerivStart[0]  -= 2.*lengthLine*xi.x()*intePs;
            l.costDerivStart[1]  -= 2.*lengthLine*xi.y()*intePs;
            l.costDerivStart[2]  -= 2.*lengthLine*xi.z()*intePs;
    
            l.costDerivEnd[0]    -= 2*lengthLine*xi.x()*intePe;
            l.costDerivEnd[1]    -= 2*lengthLine*xi.y()*intePe;
            l.costDerivEnd[2]    -= 2*lengthLine*xi.z()*intePe; 
          
        }
        
        for(unsigned i = 1; i < l.intersectionLocation.size(); i++){
            t = l.intersectionTime[i];
            PrimalPoint & PPrevious =this->primalPtList[l.intersectionLocation[i-1]];
            PrimalPoint & PNext     =this->primalPtList[l.intersectionLocation[i]];
            denominator = (PPrevious.x()-PNext.x())*(l.startPoint.x()-l.endPoint.x())+
            (PPrevious.y()-PNext.y())*(l.startPoint.y()-l.endPoint.y())+
            (PPrevious.z()-PNext.z())*(l.startPoint.z()-l.endPoint.z());
            tprime[0]=(1.-t)*(PPrevious.x()-PNext.x())/denominator;
            tprime[1]=(1.-t)*(PPrevious.y()-PNext.y())/denominator;
            tprime[2]=(1.-t)*(PPrevious.z()-PNext.z())/denominator;
            tprime[3]=(   t)*(PPrevious.x()-PNext.x())/denominator;
            tprime[4]=(   t)*(PPrevious.y()-PNext.y())/denominator;
            tprime[5]=(   t)*(PPrevious.z()-PNext.z())/denominator;
            
            double normPre;
            normPre = std::pow(PPrevious.x(),2) + std::pow(PPrevious.y(),2) + std::pow(PPrevious.z(),2);
            normPre -= std::pow(PNext.x(),2) + std::pow(PNext.y(),2) + std::pow(PNext.z(),2);
            Point p = l.getPointFromTime(t);
            double dot = p.x()*(PPrevious.x() - PNext.x()) + p.y()*(PPrevious.y() - PNext.y()) + p.z()*(PPrevious.z() - PNext.z());
            for (int k=0;k<3;k++){
                l.massDerivStart[k] -= (tprime[k]*lengthLine-t*dL[k])*(PPrevious.w()-PNext.w());
                l.massDerivEnd[k]   -= (tprime[k+3]*lengthLine+t*dL[k])*(PPrevious.w()-PNext.w());
                l.costDerivStart[k] -=  2*lengthLine*tprime[ k ]*dot;
                l.costDerivEnd[k]   -=  2*lengthLine*tprime[k+3]*dot;
                l.costDerivStart[k] += (tprime[k]*lengthLine-t*dL[k])*normPre; 
                l.costDerivEnd[k]   += (tprime[k+3]*lengthLine+t*dL[k])*normPre;
            }  
        }
        
        for (int k=0;k<3;k++){

            l.massDerivStart[k]*=l.weight;
            l.massDerivEnd  [k]*=l.weight;
            l.costDerivStart[k]*=l.weight;
            l.costDerivEnd  [k]*=l.weight;
            l.rhoDerivStart[k] =   dL[k]*(1./std::pow(polyline.totalLenght,2.));
            l.rhoDerivEnd[k]   =  -dL[k]*(1./std::pow(polyline.totalLenght,2.));
        }
    }
    
}