/*************************************************************************************************
 *
 * Modeling and animation (TNM079) 2007
 * Code base for lab assignments. Copyright:
 *   Gunnar Johansson (gunnar.johansson@itn.liu.se)
 *   Ken Museth (ken.museth@itn.liu.se)
 *   Michael Bang Nielsen (bang@daimi.au.dk)
 *   Ola Nilsson (ola.nilsson@itn.liu.se)
 *   Andreas Sderstrm (andreas.soderstrom@itn.liu.se)
 *
 *************************************************************************************************/
#include "QuadricDecimationMesh.h"

const QuadricDecimationMesh::VisualizationMode QuadricDecimationMesh::QuadricIsoSurfaces = NewVisualizationMode("Quadric Iso Surfaces");

void QuadricDecimationMesh::Initialize()
{
  // Allocate memory for the quadric array
  unsigned int numVerts = mVerts.size();
  mQuadrics.reserve(numVerts);
  std::streamsize width = std::cerr.precision(); // store stream precision
  for (unsigned int i = 0; i < numVerts; i++) {

    // Compute quadric for vertex i here
    mQuadrics.push_back(createQuadricForVert(i));


    // Calculate initial error, should be numerically close to 0

    Vector3<float> v0 = mVerts[i].pos;
    Vector4<float> v(v0[0],v0[1],v0[2],1);
    Matrix4x4<float> m = mQuadrics.back();

    float error = v*(m*v);
    //std::cerr << std::scientific << std::setprecision(2) << error << " ";
  }
  std::cerr << std::setprecision(width) << std::fixed; // reset stream precision

  // Run the initialize for the parent class to initialize the edge collapses
  DecimationMesh::Initialize();
}

/*! \lab2 Implement the computeCollapse here */
/*!
 * \param[in,out] collapse The edge collapse object to (re-)compute, DecimationMesh::EdgeCollapse
 */
void QuadricDecimationMesh::computeCollapse(EdgeCollapse * collapse)
{
  // Compute collapse->position and collapse->cost here
  // based on the quadrics at the edge endpoints

	//Matrix4x4<float> Q1 = mQuadrics.at(e(collapse->halfEdge).vert);
	//Matrix4x4<float> Q2 = mQuadrics.at(e(e(collapse->halfEdge).pair).vert);

	//Vector4<float> newPos = Vector4<float>(0.0f, 0.0f, 0.0f, 1.0f);
	//
	//Matrix4x4<float> Q12 = Q1 + Q2;
	//Matrix4x4<float> Q12copy(Q12);
	//Q12copy(3,0) = 0; Q12copy(3,1) = 0;Q12copy(3,2) = 0; Q12copy(3,3) = 1;

	//newPos = Q12copy.Inverse()*newPos;

	//collapse->position = newPos;
	//collapse->cost = newPos*(Q12*newPos);

  std::cerr << "computeCollapse in QuadricDecimationMesh is implemented, maybe\n";
}

/*! After each edge collapse the vertex properties need to be updated */
void QuadricDecimationMesh::updateVertexProperties(unsigned int ind)
{
  DecimationMesh::updateVertexProperties(ind);
  mQuadrics[ind] = createQuadricForVert(ind);
}

/*!
 * \param[in] indx vertex index, points into HalfEdgeMesh::mVerts
 */
Matrix4x4<float> QuadricDecimationMesh::createQuadricForVert(unsigned int indx) const{
  float q[4][4] = {{0,0,0,0},
                   {0,0,0,0},
                   {0,0,0,0},
                   {0,0,0,0}};
  Matrix4x4<float> Q(q);

  //// The quadric for a vertex is the sum of all the quadrics for the adjacent faces

  //std::vector<unsigned int> nFaces = FindNeighborFaces(indx);

  //for(const auto& currFace : nFaces){
	 //Q += createQuadricForFace(currFace);
  //}

  return Q;
}

/*!
 * \param[in] indx face index, points into HalfEdgeMesh::mFaces
 */
Matrix4x4<float> QuadricDecimationMesh::createQuadricForFace(unsigned int indx) const{

  // Calculate the quadric (outer product of plane parameters) for a face
  // here using the formula from Garland and Heckbert

	//Vector3<float> norm = mFaces.at(indx).normal;
	//Vector3<float> vert = v(e(mFaces.at(indx).edge).vert).pos;

	//float a = norm[0];
	//float b = norm[1];
	//float c = norm[2];
	//float d = -a*vert[0] - b*vert[1] - c*vert[2];

	////Vector4<float> p = Vector4<float>(norm[0], norm[1], norm[2], d);

	//float m[4][4];
	//m[0][0] = a*a; m[0][1] = a*b; m[0][2] = a*c; m[0][3] = a*d;
	//m[1][0] = a*b; m[1][1] = b*b; m[1][2] = b*c; m[1][3] = b*d;
	//m[2][0] = a*c; m[2][1] = c*b; m[2][2] = c*c; m[2][3] = c*d;
	//m[3][0] = a*d; m[3][1] = d*b; m[3][2] = d*c; m[3][3] = d*d;

	//Matrix4x4<float> Kp(m);

	  float q[4][4] = {{0,0,0,0},
                   {0,0,0,0},
                   {0,0,0,0},
                   {0,0,0,0}};
  Matrix4x4<float> Kp(q);

  return Kp;
}


void QuadricDecimationMesh::Render()
{
  DecimationMesh::Render();

  glEnable(GL_LIGHTING);
  glMatrixMode(GL_MODELVIEW);

  if (mVisualizationMode == QuadricIsoSurfaces)
    {
      // Apply transform
      glPushMatrix(); // Push modelview matrix onto stack

      // Implement the quadric visualization here
      std::cout << "Quadric visualization not implemented" << std::endl;

      // Restore modelview matrix
      glPopMatrix();
    }
}
