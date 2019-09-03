#include <stdexcept>
#include "polyline.h"


Line::Line(Point startPoint,Point endPoint,double weight){
    this->startPoint=startPoint;
    this->endPoint=endPoint;
    this->weight=weight;
    this->lengthLine = std::sqrt(std::pow(this->startPoint.x() - this->endPoint.x(),2) 
                                +std::pow(this->startPoint.y() - this->endPoint.y(),2) 
                                +std::pow(this->startPoint.z() - this->endPoint.z(),2));
    this->flush();
}
void Line::flush(){
    this->intersectionLocation.clear();
    this->intersectionTime.clear();
    this->massSeen =0.;
    this->cost =0.;
    this->massDerivStart = std::vector<double>(3);
    this->massDerivEnd   = std::vector<double>(3);
    this->costDerivStart = std::vector<double>(3);
    this->costDerivEnd   = std::vector<double>(3);
    this->rhoDerivStart  = std::vector<double>(3);
    this->rhoDerivEnd    = std::vector<double>(3);
    this->infoBarycenter = std::list<double>();
}

Point Line::getPointFromTime(double time){
    if(time < 0. || time > 1.){
        throw std::invalid_argument( "The time given is outside [0,1]");
    }
    return Point(this->startPoint.x()+time*(this->endPoint.x()-this->startPoint.x()),
                 this->startPoint.y()+time*(this->endPoint.y()-this->startPoint.y()),
                 this->startPoint.z()+time*(this->endPoint.z()-this->startPoint.z()));                 
}

Polyline::Polyline(double* positionX,double* positionY,double* positionZ,double* pWeight,int n){
    for(int i = 0; i < n-1 ; i++){
        this->lines.push_back(Line(Point(positionX[i],positionY[i],positionZ[i]),
                                   Point(positionX[i+1],positionY[i+1],positionZ[i+1]),
                                   pWeight[i]));
    }
    this->totalLenght = 0;
    for( auto & l : this->lines){
        totalLenght += l.lengthLine;
    }
}

Polyline::Polyline(std::list<Line>::iterator beger ,std::list<Line>::iterator ender){
    std::copy(beger, ender, std::back_inserter(this->lines));
}


Polyline::Polyline(double* positionX,double* positionY,double* positionZ,double* pWeight,int n,bool fredHDRised){
    for(int i = 0; i < n ; i++){
        this->lines.push_back(Line(Point(positionX[2*i],positionY[2*i],positionZ[2*i]),
                                   Point(positionX[2*i+1],positionY[2*i+1],positionZ[2*i+1]),
                                   pWeight[i]));
    }
    this->totalLenght = 0.0;
    for( auto & l : this->lines){
        totalLenght += l.lengthLine;
    }
}

void Polyline::flush(){
    for( auto &line : this->lines ){
        line.flush();
    }
}

void PrimalPoint::flush(){
    this->mass = 0.;
    this->cost = 0.;
    this->barycenter = std::vector<double>(3);  
    this->adjPrimalPoints.clear();
}

