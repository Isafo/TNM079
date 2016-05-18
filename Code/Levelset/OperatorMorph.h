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
#ifndef __operatormorph_h__
#define __operatormorph_h__

#include "Levelset/LevelSetOperator.h"

/*! \brief A level set operator that does morphing
 */
//! \lab5 Implement morphing
class OperatorMorph : public LevelSetOperator
{
protected :

  const Implicit * mTarget;

public :

  OperatorMorph(LevelSet * LS, const Implicit * target) : LevelSetOperator(LS), mTarget(target) { }

	virtual float ComputeTimestep()
	{
		float max = mTarget->GetValue(0,0,0);

		for(int i = 0; i < 10; i++){
			for(int j = 0; j < 10; j++){
				for(int k = 0; k < 10; k++){
					if(mTarget->GetValue(i,j,k) > max)
						max = mTarget->GetValue(i,j,k); 
				}
			}
		}
		// Compute and return a stable timestep
		//return (mLS->GetDx()/std::abs(mTarget->GetValue(x,y,z)))/2.0f;
		return mLS->GetDx() / std::abs(max);
	}

  virtual void Propagate(float time)
  {
    // Propagate level set with stable timestep dt
    // until requested time is reached
    for (float elapsed = 0; elapsed < time;) {

      // Determine timestep for stability
      float dt = ComputeTimestep();

      if (dt > time-elapsed)
        dt = time-elapsed;
      elapsed += dt;

      // Integrate level set function in time using Euler integration
      IntegrateEuler(dt);
    }
  }


virtual float Evaluate(unsigned int i, unsigned int j, unsigned int k)
{
	float x = i, y = j, z = k;

	mLS->TransformGridToWorld(x,y,z);
	float ddx2, ddy2, ddz2;
	float targetF = mTarget->GetValue(x,y,z);
	Godunov(i,j,k, targetF, ddx2, ddy2, ddz2);
	// Compute the rate of change (dphi/dt)
	return (targetF) * std::sqrt(ddx2 + ddy2 + ddz2);
}


};

#endif
