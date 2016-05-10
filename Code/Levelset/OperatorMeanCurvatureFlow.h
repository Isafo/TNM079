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
#ifndef __operatormeancurvatureflow_h__
#define __operatormeancurvatureflow_h__

#include "Levelset/LevelSetOperator.h"

/*! \brief A level set operator that does mean curvature flow.
 *
 * This class implements level set propagation in the normal direction
 * as defined by the mean curvature flow \f$\kappa\f$ in the following PDE
 *
 *  \f[
 *  \dfrac{\partial \phi}{\partial t} + \alpha \kappa|\nabla \phi| = 0
 *  \f]
 */
//! \lab4 Implement mean curvature flow
class OperatorMeanCurvatureFlow : public LevelSetOperator
{
protected:
  //! Scaling parameter, affects time step constraint
  float mAlpha;
public :

  OperatorMeanCurvatureFlow(LevelSet * LS, float alpha=.9f)
    : LevelSetOperator(LS), mAlpha(alpha) { }

  virtual float ComputeTimestep()
  {
    // Compute and return a stable timestep
	  return (std::pow(mLS->GetDx(), 2.0f)/(6*mAlpha))/2.0f;
  }

  virtual void Propagate(float time)
  {
    // Determine timestep for stability
    float dt = ComputeTimestep();

    // Propagate level set with stable timestep dt
    // until requested time is reached
    for (float elapsed = 0; elapsed < time;) {

      if (dt > time-elapsed)
        dt = time-elapsed;
      elapsed += dt;

      IntegrateEuler(dt);
      //IntegrateRungeKutta(dt);
    }
  }


  virtual float Evaluate(unsigned int i, unsigned int j, unsigned int k)
  {
    // Compute the rate of change (dphi/dt)

	  float dx = mLS->DiffXpm(i,j,k); float dx2 = dx*dx;
	  float dy = mLS->DiffYpm(i,j,k); float dy2 = dy*dy;
	  float dz = mLS->DiffZpm(i,j,k); float dz2 = dz*dz;

	  float dyz = mLS->Diff2YZpm(i,j,k);
	  float dxz = mLS->Diff2ZXpm(i,j,k);
	  float dxy = mLS->Diff2XYpm(i,j,k);

	  float dxx = mLS->Diff2Xpm(i,j,k);
	  float dyy = mLS->Diff2Ypm(i,j,k);
	  float dzz = mLS->Diff2Zpm(i,j,k);

	  float denom = (2.0f*std::pow(dx2 + dy2 + dz2, 1.5f));
	  

	  float curvature = ( (dx2*(dyy + dzz) - 2*dy*dz*dyz)/denom ) +
						( (dy2*(dxx + dzz) - 2*dx*dz*dxz)/denom ) + 
						( (dz2*(dxx + dyy) - 2*dx*dy*dxy)/denom );

	  return mAlpha*curvature*std::sqrt(dx2 + dy2 + dz2);
  }
};

#endif
