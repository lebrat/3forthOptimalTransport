#ifndef LAG3D_H
#define LAG3D_H

#include "polyline.h"
#include "optim.h"

class LaguerreCell3D{
   
    int nPoints;
    bool CoordInstanciated;
    bool PsiInstanciated;
    bool OriginInstanciated;
    double ComputeEnergy(const Point& P,const PrimalPoint& Pi);
    void SetOrigin();
    void flush();
    void GetdualInfo3d(); // Get the geometry of the Tessalation in 3d
    void GetdualInfo2d(); // Get the geometry of the Tessalation in 2d
    public:
        GenericResult FindNextIntersection(const Point& P0,const Point& P1, double actualTime, int actualLoc,int previousLoc);
        Point origin;
        int originLocation;
        Rt T;
        int Locate(const Point& P,const Point& surePoint,int sureLocation);
        int dimension(){return this->T.dimension();};
        std::vector<PrimalPoint> primalPtList;
        LaguerreCell3D();
        LaguerreCell3D(const LaguerreCell3D& other);
        void setCoordinates(double* positionX,double* positionY,double* positionZ,int nZ);
        void ComputeTesselation3D(double* weight,int nW);
        bool IsPointInCell(const Point& P,const PrimalPoint& Pi);
        void Intersect(Polyline& polyline);
        void Intersect2(Polyline& polyline);
        matCOO Intersex(Polyline& polyline);
        void IntersectParallel(Polyline& polyline,int nbThread);
        void ComputePolyLineDer(Polyline& polyline);
        int getNPoint(){return this->nPoints;};
        // void intersectionJob(LaguerreCell3D * L3,Polyline * polyline);
};



#endif