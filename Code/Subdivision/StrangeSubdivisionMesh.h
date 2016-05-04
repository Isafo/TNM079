#ifndef _strange_dubdivmesh_
#define _strange_dubdivmesh_

#include "AdaptiveLoopSubdivisionMesh.h"

class StrangeSubdivisionMesh : public AdaptiveLoopSubdivisionMesh
{
public:
  virtual void Subdivide() {
    // ....
    AdaptiveLoopSubdivisionMesh::Subdivide();
  }

protected:
  bool Subdividable(unsigned int fi){
    // Every 4th face is not subdividable - kinda strange!
    // Do something more interesting...

	   // find neighbor faces
    unsigned int f1, f2, f3;
    EdgeIterator eit = GetEdgeIterator( f(fi).edge );
    f1 = eit.Pair().GetEdgeFaceIndex(); eit.Pair();
    f2 = eit.Next().Pair().GetEdgeFaceIndex(); eit.Pair();
    f3 = eit.Next().Pair().GetEdgeFaceIndex();

	float angleSum = acos(f(fi).normal *  f(f1).normal) +
					 acos(f(fi).normal *  f(f2).normal) +
					 acos(f(fi).normal *  f(f3).normal);
	  
    return (angleSum > 0.42359877f ) ;
  }

};

#endif
