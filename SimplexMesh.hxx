// =========================================================================
// @author Leonardo Florez-Valencia (florez-l@javeriana.edu.co)
// =========================================================================
#ifndef __SimplexMesh__hxx__
#define __SimplexMesh__hxx__

#include <CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h>

// -------------------------------------------------------------------------
template< class _TKernel >
SimplexMesh< _TKernel >::
SimplexMesh( )
  : Superclass( )
{
  static const TScalar _1_3 = TScalar( 1 ) / TScalar( 3 );

  this->m_Alpha = 0.1;
  this->m_Gamma = 0.9;
  this->m_DesiredE = TVector( _1_3, _1_3, _1_3 );
  this->m_DesiredPhi = 0;
}

// -------------------------------------------------------------------------
template< class _TKernel >
SimplexMesh< _TKernel >::
SimplexMesh( const Superclass& triang )
  : Superclass( )
{
  static const TScalar _1_3 = TScalar( 1 ) / TScalar( 3 );

  this->m_Alpha = 0.1;
  this->m_Gamma = 0.9;
  this->m_DesiredE = TVector( _1_3, _1_3, _1_3 );
  this->m_DesiredPhi = 0;

  // Compute centers
  std::vector< TPoint > points( triang.number_of_faces( ) );
  for( typename Superclass::Face_index fIdx: triang.faces( ) )
  {
    THEIdx hIdx = triang.halfedge( fIdx );
    THEIdx hInit = hIdx;
    TVector bc( 0, 0, 0 );
    do
    {
      bc += triang.point( triang.source( hIdx ) ) - CGAL::ORIGIN;
      hIdx = triang.next( hIdx );
    } while( hIdx != hInit );
    bc /= 3;
    points[ int( fIdx ) ] = CGAL::ORIGIN + bc;
  } // end for

  // Assign vertices
  std::vector< TVIdx > vertices;
  for( const TPoint& p: points )
    vertices.push_back( this->add_vertex( p ) );

  // Create facets
  for( TVIdx vIdx: triang.vertices( ) )
  {
    THEIdx hIdx = triang.opposite( triang.halfedge( vIdx ) );
    THEIdx hInit = hIdx;
    std::deque< TVIdx > face_idx;
    do
    {
      face_idx.push_front( vertices[ int( triang.face( hIdx ) ) ] );
      hIdx = triang.next_around_source( hIdx );
    } while( hIdx != hInit );
    this->add_face( face_idx );
  } // end for

  // Update simplex parameters
  unsigned long nv = this->number_of_vertices( );
  this->m_N.resize( nv, TVector( 0, 0, 0 ) );
  this->m_Phi.resize( nv, 0 );
  this->m_H.resize( nv, 0 );
  this->m_E.resize( nv, TVector( 0, 0, 0 ) );
  this->m_PreviousForces.resize( nv, TVector( 0, 0, 0 ) );
  this->m_Fint.resize( nv, TVector( 0, 0, 0 ) );
  this->_ComputeSimplexParameters( );
}

// -------------------------------------------------------------------------
template< class _TKernel >
void SimplexMesh< _TKernel >::
MoveVertices( )
{
  TScalar g = TScalar( 1 ) - this->m_Gamma;

  this->_ComputeFint( );
  for( TVIdx vIdx: this->vertices( ) )
  {
    TVector d = this->m_PreviousForces[ int( vIdx ) ] * g;
    d += this->m_Fint[ int( vIdx ) ];

    this->point( vIdx ) += d;

    this->m_PreviousForces[ int( vIdx ) ] = d;
  } // end for
  this->_ComputeSimplexParameters( );
}

// -------------------------------------------------------------------------
template< class _TKernel >
void SimplexMesh< _TKernel >::
_Neighborhood(
  const TVIdx& vIdx, TPoint& p, TPoint& pn1, TPoint& pn2, TPoint& pn3
  )
{
  THEIdx hIdx = this->opposite( this->halfedge( vIdx ) );
  p = this->point( vIdx );
  pn1 = this->point( this->target( hIdx ) );
  hIdx = this->next_around_source( hIdx );
  pn2 = this->point( this->target( hIdx ) );
  hIdx = this->next_around_source( hIdx );
  pn3 = this->point( this->target( hIdx ) );
}

// -------------------------------------------------------------------------
template< class _TKernel >
typename SimplexMesh< _TKernel >::
TScalar SimplexMesh< _TKernel >::
_Project(
  TPoint& o,
  TPoint2& pp, TPoint2& ppn2, TPoint2& ppn3,
  const TVector& n,
  const TPoint& p, const TPoint& pn1, const TPoint& pn2, const TPoint& pn3
  )
{
  TVector t = pn2 - pn1;
  t /= std::sqrt( t.squared_length( ) );
  TVector b = CGAL::cross_product( t, n );
  b /= std::sqrt( b.squared_length( ) );
  ppn2 = TPoint2( t * ( pn2 - pn1 ), b * ( pn2 - pn1 ) );
  ppn3 = TPoint2( t * ( pn3 - pn1 ), b * ( pn3 - pn1 ) );
  pp   = TPoint2( t * ( p - pn1 )  , b * ( p - pn1 ) );

  // -- Circle
  TPoint2 o2 = CGAL::circumcenter( TPoint2( 0, 0 ), ppn2, ppn3 );
  o = pn1;
  o +=
    TVector(
      t[ 0 ] * o2[ 0 ] + b[ 0 ] * o2[ 1 ],
      t[ 1 ] * o2[ 0 ] + b[ 1 ] * o2[ 1 ],
      t[ 2 ] * o2[ 0 ] + b[ 2 ] * o2[ 1 ]
      );
  return(
    std::sqrt( CGAL::squared_radius( TPoint2( 0, 0 ), ppn2, ppn3 ) )
    );
}

// -------------------------------------------------------------------------
template< class _TKernel >
void SimplexMesh< _TKernel >::
_ComputeSimplexParameters( const TVIdx& vIdx )
{
  using _TBary =
    CGAL::Barycentric_coordinates::Triangle_coordinates_2< TKernel >;

  static const TScalar _0 = TScalar( 0 );

  unsigned int idx = int( vIdx );

  // -- Get point and neighbors
  TPoint p, pn1, pn2, pn3;
  this->_Neighborhood( vIdx, p, pn1, pn2, pn3 );

  // -- Mesh normal
  this->m_N[ idx ] = CGAL::normal( pn1, pn2, pn3 );

  // -- Sphere
  TPoint O = CGAL::circumcenter( p, pn1, pn2, pn3 );
  TScalar R = std::sqrt( CGAL::squared_radius( p, pn1, pn2, pn3 ) );

  // -- Projection to the ZX plane (ppn1 = < 0, 0, 0 > ) and circumcircle
  TPoint o;
  TPoint2 pp, ppn2, ppn3;
  TScalar r =
    this->_Project( o, pp, ppn2, ppn3, this->m_N[ idx ], p, pn1, pn2, pn3 );

  // -- Barycentric coordinates
  _TBary bary( TPoint2( 0, 0 ), ppn2, ppn3 );
  std::vector< TScalar > bc;
  bary( pp, std::inserter( bc, bc.end( ) ) );
  this->m_E[ idx ] = TVector( bc[ 0 ], bc[ 1 ], bc[ 2 ] );

  // -- Angle
  if( R != _0 )
  {
    TVector oO = o - O;
    TScalar cphi =
      TScalar( ( oO * this->m_N[ idx ] < _0 )? -1: 1 ) *
      std::sqrt( oO.squared_length( ) ) / R;
    TScalar sphi =
      TScalar( ( ( pn1 - p ) * this->m_N[ idx ] < _0 )? -1: 1 ) *
      r / R;
    this->m_Phi[ idx ] = std::atan2( sphi, cphi );
  }
  else
    this->m_Phi[ idx ] = _0;

  // -- Curvature
  this->m_H[ idx ] =
    ( r != _0 )? this->m_Phi[ idx ] / r: _0;
}

// -------------------------------------------------------------------------
template< class _TKernel >
void SimplexMesh< _TKernel >::
_ComputeSimplexParameters( )
{
  for( TVIdx vIdx: this->vertices( ) )
    this->_ComputeSimplexParameters( vIdx );
}

// -------------------------------------------------------------------------
template< class _TKernel >
void SimplexMesh< _TKernel >::
_ComputePhi(const TVIdx& vIdx ){
              float angle = 0;
              // -- Get point and neighbors 01
              TPoint p01, pn01, pn02, pn03;
              this->_Neighborhood( vIdx, p01, pn01, pn02, pn03 );
              static const TScalar _0 = TScalar( 0 );
              // -- Mesh normal
              //std::cout << 'c'<< p01 << pn01 << pn02 << pn03 << std::endl;
              // -- Sphere
              TPoint O = CGAL::circumcenter( p01, pn01, pn02, pn03 );
              TScalar R = std::sqrt( CGAL::squared_radius( p01, pn01, pn02, pn03 ) );

              // -- Projection to the ZX plane (ppn01 = < 0, 0, 0 > ) and circumcircle
              TPoint oo;
              TPoint2 pp01, ppn012, ppn013;
              TScalar rr =
                this->_Project( oo, pp01, ppn012, ppn013, CGAL::normal( pn01, pn02, pn03 ), p01, pn01, pn02, pn03 );

              // -- Angle
              if( R != _0 )
              {
                TVector oO = oo - O;
                TScalar cphi =
                  TScalar( ( oO * CGAL::normal( pn01, pn02, pn03 ) < _0 )? -1: 1 ) *
                  std::sqrt( oO.squared_length( ) ) / R;
                TScalar sphi =
                  TScalar( ( ( pn01 - p01 ) * CGAL::normal( pn01, pn02, pn03 ) < _0 )? -1: 1 ) *
                  rr / R;
                  angle = std::atan2( sphi, cphi );
              }
              else
                angle = 0;
              this->m_DesiredPhi = angle + this->m_DesiredPhi;
}
// -------------------------------------------------------------------------
template< class _TKernel >
void SimplexMesh< _TKernel >::
_ComputeFint( const TVIdx& vIdx )
{
  unsigned int idx = int( vIdx );

  THEIdx hIdx = this->opposite( this->halfedge( vIdx ) );
  THEIdx hIdx11 = this->next_around_source( hIdx );
  THEIdx hIdx12 = this->next_around_source( hIdx11 );

  // -- Get point and neighbors
  TPoint p, pn1, pn2, pn3;
  this->_Neighborhood( vIdx, p, pn1, pn2, pn3 );

   //desiredAnlge = 1 -> indica que toma los angulos actuales, calculados.
  int desiredAngle = 3;

  if( desiredAngle == 2){
    this->m_DesiredPhi = this->m_Phi[ idx ];
    _ComputePhi( this->target( hIdx ) );
    _ComputePhi( this->target( hIdx11 ) );
    _ComputePhi( this->target( hIdx12 ) );
  }
  std::vector< int > oldV;
  if( desiredAngle == 3 ){
    this->m_DesiredPhi = this->m_Phi[ idx ];
    oldV.push_back(idx);
    _ComputePhi( this->target( hIdx ) );
    oldV.push_back(hIdx);
    _ComputePhi( this->target( hIdx11 ) );
    oldV.push_back(hIdx11);
    _ComputePhi( this->target( hIdx12 ) );
    oldV.push_back(hIdx12);
    THEIdx hIdx22 = this->halfedge( this->target( hIdx ) );
    THEIdx hIdx223 = this->next_around_source( hIdx22 );
    if (std::find(oldV.begin(), oldV.end(),hIdx223)!=oldV.end()){

    } else {
    _ComputePhi( this->target( hIdx223 ) );
    oldV.push_back(hIdx223);
    }
    THEIdx hIdx224 = this->next_around_source( hIdx223 );
    if (std::find(oldV.begin(), oldV.end(),hIdx224)!=oldV.end()){

    } else {
    _ComputePhi( this->target( hIdx224 ) );
    oldV.push_back(hIdx224);
    }
    THEIdx hIdx225 = this->next_around_source( hIdx224 );
    if (std::find(oldV.begin(), oldV.end(),hIdx225)!=oldV.end()){

    } else {
    _ComputePhi( this->target( hIdx225 ) );
    oldV.push_back(hIdx225);
    }


    THEIdx hIdx33 = this->halfedge( this->target( hIdx11 ) );
    THEIdx hIdx333 = this->next_around_source( hIdx33 );
    if (std::find(oldV.begin(), oldV.end(),hIdx333)!=oldV.end()){

    } else {
    _ComputePhi( this->target( hIdx333 ) );
    oldV.push_back(hIdx333);
    }
    THEIdx hIdx334 = this->next_around_source( hIdx333 );
    if (std::find(oldV.begin(), oldV.end(),hIdx334)!=oldV.end()){

    } else {
    _ComputePhi( this->target( hIdx334 ) );
    oldV.push_back(hIdx334);
    }
    THEIdx hIdx335 = this->next_around_source( hIdx334 );
    if (std::find(oldV.begin(), oldV.end(),hIdx335)!=oldV.end()){

    } else {
    _ComputePhi( this->target( hIdx335 ) );
    oldV.push_back(hIdx335);
    }


    THEIdx hIdx44 = this->halfedge( this->target( hIdx12 ) );
    THEIdx hIdx444 = this->next_around_source( hIdx44 );
    if (std::find(oldV.begin(), oldV.end(),hIdx444)!=oldV.end()){

    } else {
    _ComputePhi( this->target( hIdx444 ) );
    oldV.push_back(hIdx444);
    }
    THEIdx hIdx445 = this->next_around_source( hIdx444 );
    if (std::find(oldV.begin(), oldV.end(),hIdx445)!=oldV.end()){

    } else {
    _ComputePhi( this->target( hIdx445 ) );
    oldV.push_back(hIdx445);
    }
    THEIdx hIdx446 = this->next_around_source( hIdx445 );
    if (std::find(oldV.begin(), oldV.end(),hIdx446)!=oldV.end()){

    } else {
    _ComputePhi( this->target( hIdx446 ) );
    oldV.push_back(hIdx446);
    }

  }
  if( desiredAngle == 0){
    this->m_DesiredPhi = 0;
  }
  if( desiredAngle == 1){
    this->m_DesiredPhi = this->m_Phi[ idx ];
  }

  // -- Projection to the ZX plane (ppn1 = < 0, 0, 0 > ) and circumcircle
  TPoint o;
  TPoint2 pp, ppn2, ppn3;
  TScalar r =
    this->_Project( o, pp, ppn2, ppn3, this->m_N[ idx ], p, pn1, pn2, pn3 );
  TScalar d =
    std::sqrt(
      ( pp - CGAL::circumcenter( TPoint2( 0, 0 ), ppn2, ppn3 ) ).
      squared_length( )
      );

      std::cout << this->m_DesiredPhi << std::endl;

      if(desiredAngle == 1){
        this->m_Fint[ idx ] = ( pn1 - p ) * this->m_DesiredE[ 0 ] +
        ( pn2 - p ) * this->m_DesiredE[ 1 ] +
        ( pn3 - p ) * this->m_DesiredE[ 2 ] +
        m_N[ idx ]  * Self::_L( r, d, this->m_DesiredPhi) * this->m_Alpha;
      }
      if (desiredAngle == 0) {
         this->m_Fint[ idx ] = ( pn1 - p ) * this->m_DesiredE[ 0 ] +
        ( pn2 - p ) * this->m_DesiredE[ 1 ] +
        ( pn3 - p ) * this->m_DesiredE[ 2 ] +
        m_N[ idx ]  * Self::_L( r, d, this->m_DesiredPhi) * this->m_Alpha;
      }
      if (desiredAngle == 2) {
        this->m_Fint[ idx ] = ( pn1 - p ) * this->m_DesiredE[ 0 ] +
        ( pn2 - p ) * this->m_DesiredE[ 1 ] +
        ( pn3 - p ) * this->m_DesiredE[ 2 ] +
        m_N[ idx ]  * Self::_L( r, d, this->m_DesiredPhi/4) * this->m_Alpha;
      }
      if (desiredAngle == 3) {
        this->m_Fint[ idx ] = ( pn1 - p ) * this->m_DesiredE[ 0 ] +
        ( pn2 - p ) * this->m_DesiredE[ 1 ] +
        ( pn3 - p ) * this->m_DesiredE[ 2 ] +
        m_N[ idx ]  * Self::_L( r, d, this->m_DesiredPhi/oldV.size()) * this->m_Alpha;
      }

}

// -------------------------------------------------------------------------
template< class _TKernel >
void SimplexMesh< _TKernel >::
_ComputeFint( )
{
  for( TVIdx vIdx: this->vertices( ) )
    this->_ComputeFint( vIdx );
}

// -------------------------------------------------------------------------
template< class _TKernel >
typename SimplexMesh< _TKernel >::
TScalar SimplexMesh< _TKernel >::
_L( const TScalar& r, const TScalar& d, const TScalar& p )
{
  const TScalar _pi2 = TScalar( 2 ) * std::atan( 1 );

  TScalar rd = ( r * r ) - ( d * d );
  TScalar tp = std::tan( p );
  TScalar s = TScalar( ( p < _pi2 )? 1: -1 );
  return(
    ( rd * tp ) / ( ( s * std::sqrt( ( r * r ) + ( rd * tp * tp ) ) ) + r )
    );
}

#endif // __SimplexMesh__hxx__

// eof - SimplexMesh.hxx
