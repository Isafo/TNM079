
#include "UniformCubicSplineSubdivisionCurve.h"

UniformCubicSplineSubdivisionCurve::UniformCubicSplineSubdivisionCurve(const std::vector<Vector3<float> > &joints,
                                                                       Vector3<float> lineColor,
                                                                       float lineWidth)
  : mCoefficients(joints), mControlPolygon(joints)
{
  this->mLineColor = lineColor;
  this->mLineWidth = lineWidth;
}


void UniformCubicSplineSubdivisionCurve::Subdivide()
{
	// Allocate space for new coefficients
	std::vector<Vector3<float> > newc;

	assert(mCoefficients.size() > 4 && "Need at least 5 points to subdivide");

	// Implement the subdivision scheme for a natural cubic spline here

	//newc.push_back(0.125f*(mCoefficients.at(0) + 6.0f*mCoefficients.at(0) + mCoefficients.at(1)));
	newc.push_back(mCoefficients.front());
	newc.push_back(0.125f*(4.0f*mCoefficients.at(0) + 4.0f*mCoefficients.at(1)));

	for (int i = 1; i < mCoefficients.size() - 1; ++i){
		newc.push_back(0.125f*(mCoefficients.at(i-1) + 6.0f*mCoefficients.at(i) + mCoefficients.at(i+1)));
		newc.push_back(0.125f*(4.0f*mCoefficients.at(i) + 4.0f*mCoefficients.at(i+1)));
	}

	//newc.push_back(0.125f*(mCoefficients.at(mCoefficients.size() - 2) + 6.0f*mCoefficients.back() + mCoefficients.back()));
	newc.push_back(mCoefficients.back());


	// If 'mCoefficients' had size N, how large should 'newc' be? Perform a check here!
	if(newc.size() != mCoefficients.size() + mCoefficients.size() - 1)
		std::cout << "onej" << std::endl;
	assert(true && "Incorrect number of new coefficients!");


	mCoefficients = newc;
}


void UniformCubicSplineSubdivisionCurve::Render()
{
  // Apply transform
  glPushMatrix(); // Push modelview matrix onto stack

  // Convert transform-matrix to format matching GL matrix format
  // Load transform into modelview matrix
  glMultMatrixf( mTransform.ToGLMatrix().GetArrayPtr() );

  mControlPolygon.Render();

  // save line point and color states
  glPushAttrib(GL_POINT_BIT | GL_LINE_BIT | GL_CURRENT_BIT);

  // draw segments
  glLineWidth(mLineWidth);
  glColor3fv(mLineColor.GetArrayPtr());
  glBegin(GL_LINE_STRIP);
  // just draw the spline as a series of connected linear segments
  for(unsigned int i = 0; i < mCoefficients.size(); i++){
    glVertex3fv( mCoefficients.at(i).GetArrayPtr() );
  }
  glEnd();

  // restore attribs
  glPopAttrib();

  glPopMatrix();

  GLObject::Render();
}

