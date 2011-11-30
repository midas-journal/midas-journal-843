#include <itkMesh.h>
#include <itkQuadEdgeMesh.h>

#include <iostream> 

#include "itkPointInCircleGeometricalPredicateFunctor.h"

//------------------------------------------------------------------------------
// test function
//
// main subroutine of the test
// for a given meshtype,
// test the given point against the triangle, either given as 3 points
// or as a triangle
// illustrate the code using the itk APIs
//
template< typename MeshType >
bool
test(
  const double& eps_x,
  const double& eps_y,
  const bool& ExactTest
  )
{
  // sanity check
  // check dimension, get out if <1 and issue a warning if >2

  // infer types from template
  typedef typename MeshType::PixelType  RealType;
  typedef typename MeshType::PointType  PointType;
  typedef typename MeshType::CellType   CellType; // abstract
  typedef typename MeshType::CellTraits CellTraits;

  typedef typename MeshType::CellIdentifier  CellIdentifier;
  typedef typename MeshType::PointIdentifier PointIdentifier;

  typedef typename CellType::CellAutoPointer CellAutoPointer; // abstract

  typedef itk::CellInterface< RealType, CellTraits > CellInterfaceType;
  typedef itk::TriangleCell< CellInterfaceType >     TriangleCellType;
  typedef typename TriangleCellType::PointIdIterator PointIdIterator;

  //--------------------------------
  //  direct method, just use points
  //--------------------------------

  bool result_1 = false;

  PointType mpa, mpb, mpc, mpd;
  mpa[0] = mpc[0] = mpc[1] = mpb[1] = 1.0;
  mpa[1] = mpb[0] =                   0.0;
  mpd[0] = eps_x;
  mpd[1] = eps_y;

  if( !ExactTest )
    {
    result_1 = IsInsideNotExact(
      mpa[0], mpa[1],
      mpb[0], mpb[1],
      mpc[0], mpc[1],
      mpd[0], mpd[1] );
    }
  else
    {
    result_1 = IsInside( mpa, mpb, mpc, mpd );
    }

  //----------------------------------
  // closer to real case, use triangle
  //----------------------------------

  bool result_2 = false;

  // create an empty mesh
  typedef typename MeshType::Pointer MeshPointerType;
  MeshPointerType mesh = MeshType::New();

  // add points
  mesh->SetPoint( 0, mpa );
  mesh->SetPoint( 1, mpb );
  mesh->SetPoint( 2, mpc );

  // add cell
  CellAutoPointer dummyAbstractCell;
  TriangleCellType * dummyCell = new TriangleCellType();
  dummyAbstractCell.TakeOwnership( dummyCell ); // polymorphism
  PointIdentifier dummyCellPoints[3] = {0,1,2};
  dummyAbstractCell->SetPointIds( dummyCellPoints );
  mesh->SetCell( 0, dummyAbstractCell ); // invalidate the cell

  // now this code below should be close to what the user will do
  CellIdentifier myTriangle = 0;
  result_2 =
    TestPointInTriangleInMesh< MeshType >( mesh, myTriangle, mpd, ExactTest );

  return( result_1 | result_2 );

}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// main
//
// Test for the given epsilon values, exact or not depending on user
// for a collection of different types
//
int main( int argc, char** argv )
{

  if( argc < 4 )
    {
    std::cout << "usage: <exe> epsilon_x epsilon_y IsExact";
    std::cout << std::endl;
    }

  //--------------
  // input parsing
  //--------------

  // stringPtr init
  char value = 0;
  char* ptr;
  ptr = &value;
  char** stringPtr;
  stringPtr = &ptr;

  double epsilon_x = std::strtod( argv[1], stringPtr );
  double epsilon_y = std::strtod( argv[2], stringPtr );
  std::cout << "Epsilon_x: " << epsilon_x << std::endl;
  std::cout << "Epsilon_y: " << epsilon_y << std::endl;

  bool TestExact = atoi(argv[3]);
  std::cout << "Exact Testing? " << TestExact << std::endl;

  //--------------
  // test data
  //--------------

  return(
      test< itk::Mesh< float,  2 > >( epsilon_x, epsilon_y, TestExact )
    | test< itk::Mesh< float,  3 > >( epsilon_x, epsilon_y, TestExact )
    | test< itk::Mesh< float,  4 > >( epsilon_x, epsilon_y, TestExact )
    | test< itk::Mesh< double, 2 > >( epsilon_x, epsilon_y, TestExact )
    | test< itk::Mesh< double, 3 > >( epsilon_x, epsilon_y, TestExact )
    | test< itk::Mesh< double, 4 > >( epsilon_x, epsilon_y, TestExact )
    | test< itk::QuadEdgeMesh< float,  2 > >( epsilon_x, epsilon_y, TestExact )
    | test< itk::QuadEdgeMesh< float,  3 > >( epsilon_x, epsilon_y, TestExact )
    | test< itk::QuadEdgeMesh< float,  4 > >( epsilon_x, epsilon_y, TestExact )
    | test< itk::QuadEdgeMesh< double, 2 > >( epsilon_x, epsilon_y, TestExact )
    | test< itk::QuadEdgeMesh< double, 3 > >( epsilon_x, epsilon_y, TestExact )
    | test< itk::QuadEdgeMesh< double, 4 > >( epsilon_x, epsilon_y, TestExact )
    );
}

