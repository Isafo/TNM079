/*************************************************************************************************
 *
 * Modeling and animation (TNM079) 2007
 * Code base for lab assignments. Copyright:
 *   Gunnar Johansson (gunnar.johansson@itn.liu.se)
 *   Ken Museth (ken.museth@itn.liu.se)
 *   Michael Bang Nielsen (bang@daimi.au.dk)
 *   Ola Nilsson (ola.nilsson@itn.liu.se)
 *   Andreas Sˆderstrˆm (andreas.soderstrom@itn.liu.se)
 *
 *************************************************************************************************/
#include "Geometry/HalfEdgeMesh.h"

//#include "math.h"

const unsigned int HalfEdgeMesh::BORDER = (std::numeric_limits<unsigned int>::max)();
const unsigned int HalfEdgeMesh::UNINITIALIZED = (std::numeric_limits<unsigned int>::max)()-1;

HalfEdgeMesh::HalfEdgeMesh()
{
}

HalfEdgeMesh::~HalfEdgeMesh()
{
}

/*! \lab1 Implement the addFace */
/*!
 * \param[in] v1 vertex 1, Vector3<float>
 * \param[in] v2 vertex 2, Vector3<float>
 * \param[in] v3 vertex 3, Vector3<float>
 */
bool HalfEdgeMesh::AddFace(const std::vector<Vector3<float> > &verts){
  // Add the vertices of the face/triangle
  unsigned int vI1, vI2, vI3;
  AddVertex(verts.at(0), vI1);
  AddVertex(verts.at(1), vI2);
  AddVertex(verts.at(2), vI3);

  unsigned int e1I, e1IP, e2I, e2IP, e3I, e3IP;
  // Add all half-edge pairs
  AddHalfEdgePair(vI1, vI2, e1I, e1IP);
  AddHalfEdgePair(vI2, vI3, e2I, e2IP);
  AddHalfEdgePair(vI3, vI1, e3I, e3IP);

  // Connect inner ring
  e(e1I).next = e2I;
  e(e2I).next = e3I;
  e(e3I).next = e1I;
  e(e1I).prev = e3I;
  e(e3I).prev = e2I;
  e(e2I).prev = e1I;

  // Finally, create the face, don't forget to set the normal (which should be normalized)
  Face tmpFace;
  tmpFace.edge = e1I;
  
  Vector3<float> tmpVec1 = verts.at(1) - verts.at(0);
  Vector3<float> tmpVec2 = verts.at(2) - verts.at(0);

  tmpFace.normal = Cross(tmpVec1, tmpVec2).Normalize();

  mFaces.push_back(tmpFace);

  // All half-edges share the same left face (previously added)
  unsigned int tmpInt = mFaces.size() - 1;
  e(e1I).face = tmpInt;
  e(e2I).face = tmpInt;
  e(e3I).face = tmpInt;

  // Optionally, track the (outer) boundary half-edges
  // to represent non-closed surfaces

  //e(e1IP).face = -1;
  //e(e2IP).face = -1;
  //e(e3IP).face = -1;

  return true;
}


/*!
 * \param [in] v the vertex to add, Vector3<float>
 * \param [out] indx  the index to the vertex, unsigned int
 * \return a bool indicating whether the HalfEdgeMesh::Vertex was successfully inserted (true) or already existed (false)
 */
bool HalfEdgeMesh::AddVertex(const Vector3<float> & v, unsigned int &indx){
  std::map<Vector3<float>,unsigned int>::iterator it = mUniqueVerts.find(v);
  if (it != mUniqueVerts.end()){
    indx = (*it).second; // get the index of the already existing vertex
    return false;
  }

  mUniqueVerts[v] = indx = GetNumVerts(); // op. [ ] constructs a new entry in map
  Vertex vert;
  vert.pos = v;
  mVerts.push_back(vert); // add it to the vertex list

  return true;
}

/*!
 * Inserts a half edge pair between HalfEdgeMesh::Vertex pointed to by v1 and v2. The first HalfEdgeMesh::HalfEdge (v1->v2) is the inner one, and the second (v2->v1) is the outer.
 * \param [in] v1 vertex 1, Vector3<float>
 * \param [in] v2 vertex 2, Vector3<float>
 * \param [out] indx1  the index to the half-edge from v1 to v2, unsigned int
 * \param [out] indx2  the index to the half-edge from v2 to v1, unsigned int
 * \return a bool indicating whether the half-edge pair was successfully inserted (true) or already existed (false)
 */
bool HalfEdgeMesh::AddHalfEdgePair(unsigned int v1, unsigned int v2, unsigned int &indx1, unsigned int &indx2)
{
   std::map<OrderedPair, unsigned int>::iterator it = mUniqueEdgePairs.find(OrderedPair(v1, v2));
  if (it != mUniqueEdgePairs.end()){
    indx1 = it->second;
    indx2 = e(it->second).pair;
    if(v1 != e(indx1).vert ){
      std::swap(indx1, indx2); // sort correctly
    }
    return false;
  }

  // If not found, calculate both half-edges indices
  indx1 = mEdges.size();
  indx2 = indx1+1;

  // Create edges and set pair index
  HalfEdge edge1, edge2;
  edge1.pair = indx2;
  edge2.pair = indx1;

  // Connect the edges to the verts
  edge1.vert = v1;
  edge2.vert = v2;

  // Connect the verts to the edges
  v(v1).edge = indx1;
  v(v2).edge = indx2;

  // Store the edges in mEdges
  mEdges.push_back(edge1);
  mEdges.push_back(edge2);

  // Store the first edge in the map as an OrderedPair
  OrderedPair op(v1, v2);
  mUniqueEdgePairs[op] = indx1; // op. [ ] constructs a new entry in map, ordering not important
  // sorting done when retrieving

  return true;
}

/*! \lab1 HalfEdgeMesh Implement the MergeAdjacentBoundaryEdge */
/*!
 * Merges the outer UNINITIALIZED/BORDER to an already set inner half-edge.
 * \param [in] indx the index of the INNER half-edge, unsigned int
 */
void HalfEdgeMesh::MergeOuterBoundaryEdge(unsigned int innerEdge)
{
  // Add your code here
  // 1. Merge first loop (around innerEdge->vert)
  // 2. Find leftmost edge, last edge counter clock-wise
  // 3. Test if there's anything to merge
    // 3a. If so merge the gap
    // 3b. And set border flags
  // 4. Merge second loop (around innerEdge->pair->vert)
}

/*! Proceeds to check if the mesh is valid. All indices are inspected and
 * checked to see that they are initialized. The method checks: mEdges, mFaces and mVerts.
 * Also checks to see if all verts have a neighborhood using the findNeighbourFaces method.
 */
void HalfEdgeMesh::Validate()
{
  std::vector<HalfEdge>::iterator iterEdge = mEdges.begin();
  std::vector<HalfEdge>::iterator iterEdgeEnd = mEdges.end();
  while (iterEdge != iterEdgeEnd) {
    if ((*iterEdge).face == UNINITIALIZED ||
        (*iterEdge).next == UNINITIALIZED ||
        (*iterEdge).pair == UNINITIALIZED ||
        (*iterEdge).prev == UNINITIALIZED ||
        (*iterEdge).vert == UNINITIALIZED)
        std::cerr << "HalfEdge " << iterEdge - mEdges.begin() << " not properly initialized" << std::endl;

    iterEdge++;
  }
  std::cerr << "Done with edge check (checked " << GetNumEdges() << " edges)" << std::endl;

  std::vector<Face>::iterator iterTri = mFaces.begin();
  std::vector<Face>::iterator iterTriEnd = mFaces.end();
  while (iterTri != iterTriEnd) {
    if ((*iterTri).edge == UNINITIALIZED)
        std::cerr << "Tri " << iterTri - mFaces.begin() << " not properly initialized" << std::endl;

    iterTri++;
  }
  std::cerr << "Done with face check (checked " << GetNumFaces() << " faces)" << std::endl;

  std::vector<Vertex>::iterator iterVertex = mVerts.begin();
  std::vector<Vertex>::iterator iterVertexEnd = mVerts.end();
  while (iterVertex != iterVertexEnd) {
    if ((*iterVertex).edge == UNINITIALIZED)
        std::cerr << "Vertex " << iterVertex - mVerts.begin() << " not properly initialized" << std::endl;

    iterVertex++;
  }
  std::cerr << "Done with vertex check (checked " << GetNumVerts() << " vertices)" << std::endl;

  std::cerr << "Looping through triangle neighborhood of each vertex... ";
  iterVertex = mVerts.begin();
  iterVertexEnd = mVerts.end();
  int emptyCount = 0;
  std::vector<unsigned int> problemVerts;
  while (iterVertex != iterVertexEnd) {
    std::vector<unsigned int> foundFaces = FindNeighborFaces(iterVertex - mVerts.begin());
    std::vector<unsigned int> foundVerts = FindNeighborVertices(iterVertex - mVerts.begin());
    if (foundFaces.empty() || foundVerts.empty())
      emptyCount++;
    std::set<unsigned int> uniqueFaces(foundFaces.begin(), foundFaces.end());
    std::set<unsigned int> uniqueVerts(foundVerts.begin(), foundVerts.end());
        if ( foundFaces.size() != uniqueFaces.size() ||
         foundVerts.size() != uniqueVerts.size() )
      problemVerts.push_back(iterVertex - mVerts.begin());
    iterVertex++;
  }
  std::cerr << std::endl << "Done: " << emptyCount << " isolated vertices found" << std::endl;
  if(problemVerts.size()){
    std::cerr << std::endl << "Found " << problemVerts.size() << " duplicate faces in vertices: ";
    std::copy(problemVerts.begin(), problemVerts.end(), std::ostream_iterator<unsigned int>(std::cerr, ", "));
    std::cerr << "\n";
  }
  std::cerr << std::endl << "The mesh has genus " << Genus() << ", and consists of " << Shells() << " shells.\n";
}

/*! \lab1 Implement the FindNeighborVertices */
/*! Loops over the neighborhood of a vertex and collects all the vertices sorted counter clockwise.
 * \param [in] vertexIndex  the index to vertex, unsigned int
 * \return a vector containing the indices to all the found vertices.
 */
std::vector<unsigned int> HalfEdgeMesh::FindNeighborVertices(unsigned int vertexIndex) const
{
	// Collected vertices, sorted counter clockwise!
	std::vector<unsigned int> oneRing;
	
	unsigned int firstEdge = v(vertexIndex).edge;
	unsigned int currEdge = firstEdge;

	do {
		oneRing.push_back(e(e(currEdge).pair).vert);
		currEdge = e(e(currEdge).prev).pair;
	} while(currEdge != firstEdge);

  return oneRing;
}

/*! \lab1 Implement the FindNeighborFaces */
/*! Loops over the neighborhood of a vertex and collects all the faces sorted counter clockwise.
 * \param [in] vertexIndex  the index to vertex, unsigned int
 * \return a vector containing the indices to all the found faces.
 */
std::vector<unsigned int> HalfEdgeMesh::FindNeighborFaces(unsigned int vertexIndex) const
{
	// Collected faces, sorted counter clockwise!
	std::vector<unsigned int> foundFaces;
	unsigned int firstEdge = v(vertexIndex).edge;
	unsigned int currEdge = firstEdge;

	do {
		foundFaces.push_back(e(currEdge).face);
		currEdge = e(e(currEdge).prev).pair;
	} while(currEdge != firstEdge);

	return foundFaces;
}


/*! \lab1 Implement the curvature */
// gaussian curvature
float HalfEdgeMesh::VertexCurvature(unsigned int vertexIndex) const
{

  std::vector<unsigned int> oneRing = FindNeighborVertices(vertexIndex);
  assert(oneRing.size() != 0);

  Vector3<float> aVec, bVec, cVec;

  unsigned int curr, next, prev;
  const Vector3<float> &vi = mVerts.at(vertexIndex).pos;
  float angleSum = 0;
  float area = 0;
  float cotAlpha, cotBeta;

  for(unsigned int i=0; i<oneRing.size(); i++){
    // connections
    curr = oneRing.at(i);
    if(i < oneRing.size() - 1 )
      next = oneRing.at(i+1);
    else
      next = oneRing.front();

	if(i > 0)
		prev = oneRing.at(i-1);
	else
		prev = oneRing.back();

	 // find vertices in 1-ring according to figure 5 in lab text
    // next - beta
    const Vector3<float> &nextPos = mVerts.at(next).pos;
    const Vector3<float> &vj = mVerts.at(curr).pos;
	const Vector3<float> &prevPos = mVerts.at(prev).pos; 

	//aVec = vi - nextPos;
	//bVec = vj - nextPos;
	cVec = vi - vj;
	//beta = acos((aVec.Norm() + bVec.Norm() - cVec.Norm())/(2*aVec.Length()*bVec.Length()));

	//aVec = vi - prevPos;
	//bVec = vj - prevPos;
	//alpha = acos((aVec.Norm() + bVec.Norm() - cVec.Norm())/(2*aVec.Length()*bVec.Length()));
	
	cotBeta = Cotangent(vj, nextPos, vi);

	cotAlpha = Cotangent(vi, prevPos, vj);

    // compute angle and area
    angleSum +=  acos( (vj-vi)*(nextPos-vi) / ( (vj-vi).Length()*(nextPos-vi).Length() ) );
    //area += Cross((vi-vj), (nextPos-vj)).Length()*.5;
	area += (cotAlpha + cotBeta)*cVec.Norm();
  }
  return (( 2*M_PI - angleSum ) / (area)) *0.125f;
}

// mean curvature
//float HalfEdgeMesh::VertexCurvature(unsigned int vertexIndex) const
//{
//
//  std::vector<unsigned int> oneRing = FindNeighborVertices(vertexIndex);
//  assert(oneRing.size() != 0);
//
//  Vector3<float> aVec, bVec, cVec;
//
//  unsigned int curr, next, prev;
//  const Vector3<float> &vi = mVerts.at(vertexIndex).pos;
//  Vector3<float> twoGrad = Vector3<float>(0.0f, 0.0f, 0.0f);
//  float area = 0;
//  float cotAlpha, cotBeta;
//
//  for(unsigned int i=0; i<oneRing.size(); i++){
//    // connections
//    curr = oneRing.at(i);
//    if(i < oneRing.size() - 1 )
//      next = oneRing.at(i+1);
//    else
//      next = oneRing.front();
//
//	if(i > 0)
//		prev = oneRing.at(i-1);
//	else
//		prev = oneRing.back();
//
//	 // find vertices in 1-ring according to figure 5 in lab text
//    // next - beta
//    const Vector3<float> &nextPos = mVerts.at(next).pos;
//    const Vector3<float> &vj = mVerts.at(curr).pos;
//	const Vector3<float> &prevPos = mVerts.at(prev).pos; 
//
//	//aVec = vi - nextPos;
//	//bVec = vj - nextPos;
//	cVec = vi - vj;
//	//beta = acos((aVec.Norm() + bVec.Norm() - cVec.Norm())/(2*aVec.Length()*bVec.Length()));
//
//	//aVec = vi - prevPos;
//	//bVec = vj - prevPos;
//	//alpha = acos((aVec.Norm() + bVec.Norm() - cVec.Norm())/(2*aVec.Length()*bVec.Length()));
//	
//	cotBeta = Cotangent(vj, nextPos, vi);
//
//	cotAlpha = Cotangent(vi, prevPos, vj);
//
//    // compute angle and area
//	twoGrad += (cotAlpha + cotBeta)*cVec;
//    //area += Cross((vi-vj), (nextPos-vj)).Length()*.5;
//	area += (cotAlpha + cotBeta)*cVec.Norm();
//  }
//  return twoGrad.Length()/(0.5*area);
//}

float HalfEdgeMesh::FaceCurvature(unsigned int faceIndex) const
{
  // NB Assumes vertex curvature already computed
  unsigned int indx = f(faceIndex).edge;
  const EdgeIterator it = GetEdgeIterator(indx);

  const Vertex& v1 = v(it.GetEdgeVertexIndex());
  const Vertex &v2 = v(it.Next().GetEdgeVertexIndex());
  const Vertex &v3 = v(it.Next().GetEdgeVertexIndex());

  return (v1.curvature + v2.curvature + v3.curvature) / 3.f;
}

Vector3<float> HalfEdgeMesh::FaceNormal(unsigned int faceIndex) const
{
  unsigned int indx = f(faceIndex).edge;
  const EdgeIterator it = GetEdgeIterator(indx);

  const Vector3<float> &p1 = v(it.GetEdgeVertexIndex()).pos;
  const Vector3<float> &p2 = v(it.Next().GetEdgeVertexIndex()).pos;
  const Vector3<float> &p3 = v(it.Next().GetEdgeVertexIndex()).pos;

  const Vector3<float> e1 = p2-p1;
  const Vector3<float> e2 = p3-p1;
  return Cross(e1, e2).Normalize();
}

Vector3<float> HalfEdgeMesh::VertexNormal(unsigned int vertexIndex) const
{
	Vector3<float> n(0,0,0);

	std::vector<unsigned int> faces = FindNeighborFaces(vertexIndex);

	for(int i = 0; i < faces.size(); i++) {
		n += f(faces.at(i)).normal;
	}

	n.Normalize();

	return n;
}


void HalfEdgeMesh::Initialize() {
  Validate();
  Update();
}


void HalfEdgeMesh::Update() {
  // Calculate and store all differentials and area

  // First update all face normals and triangle areas
  for(unsigned int i = 0; i < GetNumFaces(); i++){
    f(i).normal = FaceNormal(i);
  }
  // Then update all vertex normals and curvature
  for(unsigned int i = 0; i < GetNumVerts(); i++){
    // Vertex normals are just weighted averages
    mVerts.at(i).normal = VertexNormal(i);
  }

  // Then update vertex curvature
  for(unsigned int i = 0; i < GetNumVerts(); i++){
    mVerts.at(i).curvature = VertexCurvature(i);
    //    std::cerr <<   mVerts.at(i).curvature << "\n";
  }

  // Finally update face curvature
  for(unsigned int i = 0; i < GetNumFaces(); i++){
    f(i).curvature = FaceCurvature(i);
  }

  std::cerr << "Area: " << Area() << ".\n";
  std::cerr << "Volume: " << Volume() << ".\n";

  // Update vertex and face colors
  if (mVisualizationMode == CurvatureVertex) {
    std::vector<Vertex>::iterator iter = mVerts.begin();
    std::vector<Vertex>::iterator iend = mVerts.end();
    float minCurvature = (std::numeric_limits<float>::max)();
    float maxCurvature = -(std::numeric_limits<float>::max)();
    while (iter != iend) {
      if (minCurvature > (*iter).curvature)  minCurvature = (*iter).curvature;
      if (maxCurvature < (*iter).curvature)  maxCurvature = (*iter).curvature;
      iter++;
    }
    std::cerr << "Mapping color based on vertex curvature with range [" << minCurvature << "," << maxCurvature << "]" << std::endl;
    iter = mVerts.begin();
    while (iter != iend) {
      (*iter).color = mColorMap->Map((*iter).curvature, minCurvature, maxCurvature);
      iter++;
    }
  }
  else if (mVisualizationMode == CurvatureFace) {
    std::vector<Face>::iterator iter = mFaces.begin();
    std::vector<Face>::iterator iend = mFaces.end();
    float minCurvature = (std::numeric_limits<float>::max)();
    float maxCurvature = -(std::numeric_limits<float>::max)();
    while (iter != iend) {
      if (minCurvature > (*iter).curvature)  minCurvature = (*iter).curvature;
      if (maxCurvature < (*iter).curvature)  maxCurvature = (*iter).curvature;
      iter++;
    }
    std::cerr << "Mapping color based on face curvature with range [" << minCurvature << "," << maxCurvature << "]" << std::endl;
    iter = mFaces.begin();
    while (iter != iend) {
      (*iter).color = mColorMap->Map((*iter).curvature, minCurvature, maxCurvature);
      iter++;
    }
  }

}


/*! \lab1 Implement the area */
float HalfEdgeMesh::Area() const
{
	float area = 0;

	unsigned int v1, v2, v3;
	Vector3<float> tempVec1, tempVec2;

	for (int i = 0; i < mFaces.size(); i++){
		v1 = e(f(i).edge).vert;
		v2 = e(e(f(i).edge).next).vert;
		v3 = e(e(f(i).edge).prev).vert;

		tempVec1 = v(v2).pos - v(v1).pos;
		tempVec2 = v(v3).pos - v(v1).pos;

		area += Cross(tempVec1, tempVec2).Length() * 0.5;

	}
	return area;
}

/*! \lab1 Implement the volume */
float HalfEdgeMesh::Volume() const
{
	float volume = 0;

  
	unsigned int v1, v2, v3;
	Vector3<float> tempVec1, tempVec2, centroid;
	float area;


  	for (int i = 0; i < mFaces.size(); i++){
		v1 = e(f(i).edge).vert;
		v2 = e(e(f(i).edge).next).vert;
		v3 = e(e(f(i).edge).prev).vert;

		centroid =  (v(v1).pos +  v(v2).pos +  v(v3).pos)/3;

		tempVec1 = v(v2).pos - v(v1).pos;
		tempVec2 = v(v3).pos - v(v1).pos;

		area = Cross(tempVec1, tempVec2).Length() * 0.5;

		volume += (centroid * f(i).normal) * area;
	}

  return volume/3.0f;
}

/*! \lab1 Calculate the number of shells  */
int HalfEdgeMesh::Shells() const
{
	std::set<unsigned int> vertexQueueSet;
	std::set<unsigned int> vertexTaggedSet;

	vertexQueueSet.insert(0);

	unsigned int currVert;
	std::vector<unsigned int> neighbourVerts;

	int shells = 0;

	while(true) { 
		do {
			currVert = (*vertexQueueSet.begin());
			vertexQueueSet.erase(vertexQueueSet.begin());

			vertexTaggedSet.insert(currVert);
		
			neighbourVerts = FindNeighborVertices(currVert);

			for(const auto& nVert : neighbourVerts) {
				if(vertexTaggedSet.find(nVert) != vertexTaggedSet.end()) 
					continue;

				vertexQueueSet.insert(nVert);
			}
		} while(!vertexQueueSet.empty());

		++shells;

		if(vertexTaggedSet.size() == mUniqueVerts.size())
			break;

		for(const auto& vert : mUniqueVerts) {
			if(vertexTaggedSet.find(vert.second) == vertexTaggedSet.end()) {
				vertexQueueSet.insert(vert.second);
				break;
			}
		}	
	}


	return shells;
}

/*! \lab1 Implement the genus */
int HalfEdgeMesh::Genus() const
{
	return (-mUniqueVerts.size() + mUniqueEdgePairs.size() - mFaces.size() + 2 * Shells() ) / 2;
}

void HalfEdgeMesh::Dilate(float amount)
{
  std::vector<Vertex>::iterator iter = mVerts.begin();
  std::vector<Vertex>::iterator iend = mVerts.end();
  while (iter != iend) {
    (*iter).pos += amount*(*iter).normal;
    iter++;
  }

  Initialize();
  Update();
}

void HalfEdgeMesh::Erode(float amount)
{
  std::vector<Vertex>::iterator iter = mVerts.begin();
  std::vector<Vertex>::iterator iend = mVerts.end();
  while (iter != iend) {
    (*iter).pos -= amount*(*iter).normal;
    iter++;
  }

  Initialize();
  Update();
}

void HalfEdgeMesh::Smooth(float amount)
{
  std::vector<Vertex>::iterator iter = mVerts.begin();
  std::vector<Vertex>::iterator iend = mVerts.end();
  while (iter != iend) {
    (*iter).pos -= amount*(*iter).normal * (*iter).curvature;
    iter++;
  }

  Initialize();
  Update();
}

void HalfEdgeMesh::Render()
{
  glEnable(GL_LIGHTING);
  glMatrixMode(GL_MODELVIEW);

  // Apply transform
  glPushMatrix(); // Push modelview matrix onto stack

  // Convert transform-matrix to format matching GL matrix format
  // Load transform into modelview matrix
  glMultMatrixf(  mTransform.ToGLMatrix().GetArrayPtr() );

  if (mWireframe)
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

  // Draw geometry
  glBegin(GL_TRIANGLES);
  const int numTriangles = GetNumFaces();
  for (int i = 0; i < numTriangles; i++){

    Face & face = f(i);

    HalfEdge* edge = &e(face.edge);

    Vertex & v1 = v(edge->vert);
    edge = &e(edge->next);

    Vertex & v2 = v(edge->vert);
    edge = &e(edge->next);

    Vertex & v3 = v(edge->vert);

    if (mVisualizationMode == CurvatureVertex) {
      glColor3fv(v1.color.GetArrayPtr());
      glNormal3fv(v1.normal.GetArrayPtr());
      glVertex3fv(v1.pos.GetArrayPtr());

      glColor3fv(v2.color.GetArrayPtr());
      glNormal3fv(v2.normal.GetArrayPtr());
      glVertex3fv(v2.pos.GetArrayPtr());

      glColor3fv(v3.color.GetArrayPtr());
      glNormal3fv(v3.normal.GetArrayPtr());
      glVertex3fv(v3.pos.GetArrayPtr());
    }
    else {
      glColor3fv(face.color.GetArrayPtr());
      glNormal3fv(face.normal.GetArrayPtr());

      glVertex3fv(v1.pos.GetArrayPtr());
      glVertex3fv(v2.pos.GetArrayPtr());
      glVertex3fv(v3.pos.GetArrayPtr());
    }

  }
  glEnd();

  // Mesh normals by courtesy of Richard Khoury
  if (mShowNormals)
  {
    glDisable(GL_LIGHTING);
    glBegin(GL_LINES);
    const int numTriangles = GetNumFaces();
    for (int i = 0; i < numTriangles; i++){

      Face & face = f(i);

      HalfEdge* edge = &e(face.edge);

      Vertex & v1 = v(edge->vert);
      edge = &e(edge->next);

      Vertex & v2 = v(edge->vert);
      edge = &e(edge->next);

      Vertex & v3 = v(edge->vert);

      Vector3<float> faceStart = (v1.pos + v2.pos + v3.pos) / 3.0f;
      Vector3<float> faceEnd = faceStart + face.normal*0.1;

      glColor3f(1,0,0); // Red for face normal
      glVertex3fv(faceStart.GetArrayPtr());
      glVertex3fv(faceEnd.GetArrayPtr());

      glColor3f(0,1,0); // Vertex normals in Green
      glVertex3fv(v1.pos.GetArrayPtr());
      glVertex3fv((v1.pos + v1.normal*0.1).GetArrayPtr());
      glVertex3fv(v2.pos.GetArrayPtr());
      glVertex3fv((v2.pos + v2.normal*0.1).GetArrayPtr());
      glVertex3fv(v3.pos.GetArrayPtr());
      glVertex3fv((v3.pos + v3.normal*0.1).GetArrayPtr());
    }
    glEnd();
    glEnable(GL_LIGHTING);

  }


  if (mWireframe)
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  // Restore modelview matrix
  glPopMatrix();

  GLObject::Render();

}

