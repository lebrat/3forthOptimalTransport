#ifndef POLYLINE_H
#define POLYLINE_H


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/config.h>
#include <limits>
#include <vector>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT                                               Weight;
typedef K::Point_3                                          Point;
typedef K::Weighted_point_3                                 Weighted_point;
typedef CGAL::Regular_triangulation_vertex_base_3<K>        Vb0;
typedef CGAL::Triangulation_vertex_base_with_info_3<int, K, Vb0> Vb;
typedef CGAL::Regular_triangulation_cell_base_3<K>          Cb;
typedef CGAL::Triangulation_data_structure_3<Vb,Cb>         Tds;
typedef CGAL::Regular_triangulation_3<K, Tds>               Rt;
typedef Rt::Vertex_handle                                   Vertex_handle;
typedef Rt::Finite_facets_iterator                          Finite_facets_iterator;
typedef Rt::All_facets_iterator                             All_facets_iterator;
typedef Rt::Cell_iterator                                   Cell_iterator;

struct GenericResult{
    double time;
    int location; 
};


class PrimalPoint{
    public:
        Point point; // Localization of the primal point
        double psi; //Value of psi for the primal point
        // OPTIMAL TRANSPORT INFO
        double mass; //Mass of the polyline in the primal cell
        double cost; //Cost of the transport in the primal cell
        std::vector<double> barycenter; //Barycentre of the polyline in the primal cell
        std::set<int> adjPrimalPoints; // List of neighbours primal points
        void flush(); //flushes all the above optimal transport info
        // GEOMETRY INFO
        int id;
        double x(void) const {return this->point.x();};
        double y(void) const {return this->point.y();};
        double z(void) const {return this->point.z();};
        double w(void) const {return this->psi;};
        void w(double w) {this->psi=w;};

        PrimalPoint(){};
        PrimalPoint(Point p,int i){
            this->point=p;
            this->id=i;
            this->flush();
        };
};

struct Line {
    Point getPointFromTime(double time);

    double computeSquaredDistance(double t, const PrimalPoint & P);
    // Geometric info
    Point startPoint;
    Point endPoint;
    double weight;
    double lengthLine;
    // OPTIMAL TRANSPORT INFO
    std::vector<int> intersectionLocation;
    std::vector<double> intersectionTime;
    double massSeen; // Total mass of points that is send to the polyline
    double cost; // Cost of the transport of the polyline
    std::vector<double> massDerivStart,massDerivEnd, costDerivStart,costDerivEnd;
    std::vector<double> rhoDerivStart,rhoDerivEnd;
    std::list<double> infoBarycenter;
    void flush(); //flushes all the above optimal transport info
    Line(){};
    Line(Point startPoint,Point endPoint,double weight);
};

struct Polyline{
    std::list<Line> lines;
    double totalLenght;
    Polyline(){};
    Polyline(double* positionX,double* positionY,double* positionZ,double* weight,int n);
    Polyline(std::list<Line>::iterator beger ,std::list<Line>::iterator ender);
    Polyline(double* positionX,double* positionY,double* positionZ,double* weight,int n,bool fredHDRised);
    void flush();  
};


#endif