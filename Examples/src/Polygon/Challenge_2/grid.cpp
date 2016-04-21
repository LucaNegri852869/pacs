#include "Polygon.hpp"
#include "grid.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

namespace Geometry
{
void Grid::read_data(std::ifstream & infile){
  unsigned int n_points(0);
  unsigned int n_poly(0);
  unsigned int PointID(0);
  unsigned int type(0);
  unsigned int PolyID(0);
  double x(0.0);
  double y(0.0);
  Vertices vertexes;
  infile>>n_points;
  infile>>n_poly;
  //Read all the coordinates of the points
  for(unsigned int i=0;i<n_points;i++){
    infile>>PointID;
    infile>>x;
    infile>>y;
    Mesh_Point.push_back(Point2D(x,y));
    }
   //Read all the polygons 
   for(unsigned int j=0;j<n_poly;j++){
     infile>>PolyID;
     infile>>type;
     switch(type){
     case(0):
     {
     	for(int k=0;k<3;k++){
     	  infile>>PointID;
     	  vertexes.push_back(Mesh_Point[PointID]);
     	  }
     	  std::shared_ptr<Triangle> tri(new Triangle(vertexes));
     	  Mesh_Poly.push_back(tri);
     	  vertexes.clear();
     	  break;
     }
     case(1):
     {
     	for(int k=0;k<4;k++){
     	  infile>>PointID;
     	  vertexes.push_back(Mesh_Point[PointID]);
     	  }
     	  std::shared_ptr<Square> sq(new Square(vertexes));
     	  Mesh_Poly.push_back(sq);
     	  vertexes.clear();
     	  break;
     }
     case(2):
     {
     	for(int k=0;k<4;k++){
     	  infile>>PointID;
     	  vertexes.push_back(Mesh_Point[PointID]);
     	  }
     	  std::shared_ptr<Polygon> abs(new Polygon(vertexes));
     	  Mesh_Poly.push_back(abs);
     	  vertexes.clear();
     	  break;
     }
     }
   }
 }
 
 
double Grid::total_area() {
  double sum(0.0);
  for(unsigned int i=0;i<Mesh_Poly.size();i++){
    sum+=Mesh_Poly[i]->area();
    }
  return sum;
}
   
     
}


