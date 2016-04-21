#ifndef HH_GRID_HH
#define HH_GRID_HH
#include <iostream>
#include <vector>
#include <array>
#include <utility>
#include <memory>
#include "Polygon.hpp"

namespace Geometry
{
class Grid{
   public:
   //default constructor
   Grid()=default;
   //default copy constructor
   Grid(Grid const &)=default;
   //assignement operator
   Grid & operator=(Grid const &)=default;
   
   //Method that reads the file
   void read_data(std::ifstream &infile);
   
   //Method that computes the sum of the areas of the polygons
   double total_area() ;
   
   private:
   //vector of shared pointers to AbsPoly
   std::vector<std::shared_ptr <AbstractPolygon> > Mesh_Poly;
   //vector of Point2D
   Vertices Mesh_Point;
   };
  
  }

#endif
