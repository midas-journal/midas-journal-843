#include <itkMatrix.h>
#include <itkTriangleCell.h>
#include <itkCellInterface.h>

#include "vnl/vnl_det.h" 
#include "vnl/vnl_matrix_fixed.h" 

//------------------------------------------------------------------------------
// Skewchuck code
//
extern "C"
  {
  double incircle(double* pa, double* pb, double* pc, double* pd);
  double orient2d(double* pa, double* pb, double* pc);
  }
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// The test functor to map skewchuck code to ITK's API
//
template< typename PointType >
bool
IsInside(
  const PointType& TrianglePoint1,
  const PointType& TrianglePoint2,
  const PointType& TrianglePoint3,
  const PointType& PointToTest )
{

  double * pa = new double[2];
  double * pb = new double[2];
  double * pc = new double[2];
  double * pd = new double[2]; 
  pa[0] = TrianglePoint1[0];
  pa[1] = TrianglePoint1[1];
  pb[0] = TrianglePoint2[0];
  pb[1] = TrianglePoint2[1];
  pc[0] = TrianglePoint3[0];
  pc[1] = TrianglePoint3[1];
  pd[0] = PointToTest[0];
  pd[1] = PointToTest[1];

  // orientation test
  double orientation = orient2d( pa, pb, pc );

  // incircle test - the result is multiplied by the orientation test result
  double det = incircle( pa, pb, pc, pd ) * orientation;
  
  // NOTE ALEX: replace by debug macro
  std::cout << "Det(M): " << det << std::endl;
 
  // NOTE ALEX: zero, which means the point is ON the circle is considered IN.
  return( det>=0?0:1 );
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
//  Non exact, ITK native version
//
template< typename RealType >
bool
IsInsideNotExact(
  const RealType& ax,
  const RealType& ay,
  const RealType& bx,
  const RealType& by,
  const RealType& cx,
  const RealType& cy,
  const RealType& dx,
  const RealType& dy )
{
  double det;

  // orientation test - determination of the sign of ad-bc
  double a = ax-cx;
  double b = ay-cy;
  double c = bx-cx;
  double d = by-cy;
  double orientation = a*d-b*c;

  typedef itk::Matrix< double, 4,4 > MatrixType;
  MatrixType M;

  M(0,0) = ax;
  M(0,1) = ay;
  M(0,2) = ax*ax+ay*ay;
  M(0,3) = 1.0;
  M(1,0) = bx;
  M(1,1) = by;
  M(1,2) = bx*bx+by*by;
  M(1,3) = 1.0;
  M(2,0) = cx;
  M(2,1) = cy;
  M(2,2) = cx*cx+cy*cy;
  M(2,3) = 1.0;
  M(3,0) = dx;
  M(3,1) = dy;
  M(3,2) = dx*dx+dy*dy;
  M(3,3) = 1.0;
 
  std::cout << "M: " << std::endl;
  std::cout << M << std::endl;

  // determinant computation - the result is multiplied by 'orientation'
  det = vnl_det( M.GetVnlMatrix() ) * orientation;
  std::cout << "Det(M): " << det << std::endl;

  return( det>=0?0:1 );
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// Convenience Functor
//
// fetch triangle from the container,
// fetch pointsIds from triangleCell
// fetch points from point container
// and test
//
// NOTE ALEX: this function still use the non exact, naive implementation
// of the predicate that should not remain in the final itk version.
// to be cleaned during the review process.
//
template<
  typename MeshType, 
  typename MeshPointerType, 
  typename CellIdentifier, 
  typename PointType >
bool
TestPointInTriangleInMesh(
  MeshPointerType myMesh,
  const CellIdentifier& TriangleId,
  const PointType& PointToTest,
  bool ExactTest
  )
{
  typedef typename MeshType::CellType   CellType; // abstract
  typedef typename CellType::CellAutoPointer CellAutoPointer; // abstract

  typedef typename MeshType::PixelType  RealType;
  typedef typename MeshType::CellTraits CellTraits;

  typedef itk::CellInterface< RealType, CellTraits > CellInterfaceType;
  typedef itk::TriangleCell< CellInterfaceType >     TriangleCellType;
  typedef typename TriangleCellType::PointIdIterator PointIdIterator;

  PointType mpa, mpb, mpc;

  CellAutoPointer cellIterator;
  if( myMesh->GetCell( TriangleId, cellIterator ) )
    {
    if( cellIterator->GetType() == 2 )
      { 
      // we have a triangle, let s get the points
      PointIdIterator pointIdIterator = cellIterator->PointIdsBegin();
      // NOTE ALEX: should check the return value
      myMesh->GetPoint( *pointIdIterator, &mpa );
      pointIdIterator++;
      myMesh->GetPoint( *pointIdIterator, &mpb );
      pointIdIterator++;
      myMesh->GetPoint( *pointIdIterator, &mpc );
      if( !ExactTest )
        {
        return IsInsideNotExact( //< RealType >(
          mpa[0], mpa[1],
          mpb[0], mpb[1],
          mpc[0], mpc[1],
          PointToTest[0], PointToTest[1] );
        }
      else
        {
        return IsInside( mpa, mpb, mpc, PointToTest);
        }
      }
    else
      {
      // oops - this is not a triangle
      return( false );
      }
    }
  else
   {
   // oops, the cell ID was not found in the container
   return( false );
   }
}
//------------------------------------------------------------------------------

