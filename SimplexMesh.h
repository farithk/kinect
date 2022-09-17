// =========================================================================
// @author Leonardo Florez-Valencia (florez-l@javeriana.edu.co)
// =========================================================================
#ifndef __SimplexMesh__h__
#define __SimplexMesh__h__

#include <vector>
#include <CGAL/Surface_mesh.h>

/**
 */
template< class _TKernel >
class SimplexMesh
  : public CGAL::Surface_mesh< typename _TKernel::Point_3 >
{
public:
  using Self       = SimplexMesh;
  using Superclass = CGAL::Surface_mesh< typename _TKernel::Point_3 >;

  using TKernel  = _TKernel;
  using TScalar  = typename TKernel::FT;
  using TPoint   = typename TKernel::Point_3;
  using TVector  = typename TKernel::Vector_3;
  using TPoint2  = typename TKernel::Point_2;
  using TVector2 = typename TKernel::Vector_2;

  using THEIdx = typename Superclass::Halfedge_index;
  using TVIdx  = typename Superclass::Vertex_index;

public:
  SimplexMesh( );
  SimplexMesh( const Superclass& triang );
  virtual ~SimplexMesh( ) = default;

  void MoveVertices( );

protected:
  void _Neighborhood(
    const TVIdx& vIdx, TPoint& p, TPoint& pn1, TPoint& pn2, TPoint& pn3
    );
  TScalar _Project(
    TPoint& o,
    TPoint2& pp, TPoint2& ppn2, TPoint2& ppn3,
    const TVector& n,
    const TPoint& p, const TPoint& pn1, const TPoint& pn2, const TPoint& pn3
    );

  void _ComputeSimplexParameters( const TVIdx& vIdx );
  void _ComputeSimplexParameters( );

  void _ComputeFint( const TVIdx& vIdx );
  void _ComputeFint( );

  void _ComputePhi( const TVIdx& vIdx );

  static TScalar _L( const TScalar& r, const TScalar& d, const TScalar& p );

protected:
  std::vector< TVector > m_N;
  std::vector< TScalar > m_Phi;
  std::vector< TScalar > m_H;
  std::vector< TVector > m_E;

  TScalar m_Alpha;
  TScalar m_Gamma;
  TVector m_DesiredE;
  TScalar m_DesiredPhi;
  std::vector< TVector > m_PreviousForces;
  std::vector< TVector > m_Fint;
};

#include "SimplexMesh.hxx"

#endif // __SimplexMesh__h__

// eof - SimplexMesh.h
