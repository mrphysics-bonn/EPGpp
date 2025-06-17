/** @file epg.cpp
 *  @brief Implementation of Extended Phase Graph Signal Simulations
 *
 * Author: tony.stoecker@dzne.de
 * Date  : August 20, 2013
 */

#include "epg.h"

/****************************************************************/
EPG::EPG (const double &M0 , const double &T1 , const double &T2, const double &TR ) {
	SetParameters(M0,T1,T2,TR); 
	m_verbose = false;
};

/****************************************************************/
EPG::~EPG (        ) {	DeleteStates(); };

/****************************************************************/
EPG::EPG (const EPG &epg) { Copy(epg); };

/****************************************************************/
void EPG::Copy (const EPG &epg){
  
	SetParameters(epg.GetM0(),epg.GetT1(),epg.GetT2(),epg.GetTR()); 
	m_step    = epg.GetStep();
	m_verbose = epg.GetVerbose();
	
	DeleteStates();
	int i;
	for (i=0; i<epg.m_FaRe.size(); ++i)  { m_FaRe.push_back(epg.m_FaRe[i]); m_FaIm.push_back(epg.m_FaIm[i]); }
	for (i=0; i<epg.m_FbRe.size(); ++i)  { m_FbRe.push_back(epg.m_FbRe[i]); m_FbIm.push_back(epg.m_FbIm[i]); }
	for (i=0; i<epg.m_ZaRe.size(); ++i)  { m_ZaRe.push_back(epg.m_ZaRe[i]); m_ZaIm.push_back(epg.m_ZaIm[i]); }
	for (i=0; i<epg.m_ZbRe.size(); ++i)  { m_ZbRe.push_back(epg.m_ZbRe[i]); m_ZbIm.push_back(epg.m_ZbIm[i]); }

};

/****************************************************************/
EPG& EPG::operator= (const EPG &epg){
  Copy(epg);
  return *this;
};

/****************************************************************/
bool EPG::operator== (const EPG &epg){

	if ( m_M0 != epg.m_M0 ) return false; 
	if ( m_T1 != epg.m_T1 ) return false; 
	if ( m_T2 != epg.m_T2 ) return false; 
	
	if ( m_step != epg.m_step ) return false; 
	
	if ( m_verbose != epg.m_verbose ) return false; 
	
	if ( m_FaRe != epg.m_FaRe ) return false; 
	if ( m_FbRe != epg.m_FbRe ) return false; 
	if ( m_ZaRe != epg.m_ZaRe ) return false; 
	if ( m_ZbRe != epg.m_ZbRe ) return false; 

	if ( m_FaIm != epg.m_FaIm ) return false; 
	if ( m_FbIm != epg.m_FbIm ) return false; 
	if ( m_ZaIm != epg.m_ZaIm ) return false; 
	if ( m_ZbIm != epg.m_ZbIm ) return false; 
	
	return true;
};

/****************************************************************/
bool EPG::operator!= (const EPG &epg){
	return !(*this == epg);
};

/****************************************************************/
void EPG::SetParameters (const double &M0, const double &T1, const double &T2, const double &TR ) {
	m_M0 = M0;
	m_T1 = T1;
	m_T2 = T2;
	m_TR = TR;
	m_E1 = exp(-TR/m_T1);
	m_E2 = exp(-TR/m_T2);
	Equilibrium();
};

/****************************************************************/
void EPG::DeleteStates(){
	m_FaRe.clear();
	m_FaIm.clear();
	m_FbRe.clear();
	m_FbIm.clear();
	m_ZaRe.clear();
	m_ZaIm.clear();
	m_ZbRe.clear();
	m_ZbIm.clear();
 };

/****************************************************************/
void EPG::Equilibrium(){
	DeleteStates();
	m_FaRe.push_back(0.0);
	m_FaIm.push_back(0.0);
	m_FbRe.push_back(0.0);
	m_FbIm.push_back(0.0);
	m_ZaRe.push_back(0.0);
	m_ZaIm.push_back(0.0);
	m_ZbRe.push_back(0.0);
	m_ZbIm.push_back(0.0);
	m_ZbRe[0]=m_M0;
  	m_step=0;
};

/****************************************************************/
void EPG::NullTransverse(){
	for (int i=0; i<m_FaRe.size(); ++i) { m_FaRe[i]=0.0; m_FaIm[i]=0.0;  m_FbRe[i]=0.0; m_FbIm[i]=0.0; }
};

/****************************************************************/
void EPG::SetStep(int val) {
//
	int SFa = m_FaRe.size();
	int SFb = m_FbRe.size();
	int SZa = m_ZaRe.size();
	int SZb = m_ZbRe.size();
	if (val>SZa) return;
// 	
	m_step=val; 
	int i;

	for ( i=2*m_step+1; i<SFa; ++i ) { m_FaRe[i]=0.0; m_FaIm[i]=0.0;}
	for ( i=2*m_step+1; i<SFb; ++i ) { m_FbRe[i]=0.0; m_FbIm[i]=0.0;}
	for ( i=  m_step  ; i<SZa; ++i ) { m_ZaRe[i]=0.0; m_ZaIm[i]=0.0;}
	for ( i=  m_step  ; i<SZb; ++i ) { m_ZbRe[i]=0.0; m_ZbIm[i]=0.0;}
};

/****************************************************************/
void   EPG::GetFaState ( double* real, double* imag, const int &num ) const {
	int N = m_step-num-1;
	if ( N >= 0 && N < 2*m_step+1 ) { *real = m_FaRe[N]; *imag = m_FaIm[N]; }
};

void   EPG::GetFbState ( double* real, double* imag, const int &num ) const {
	int N = m_step-num-1;
	if ( N >= 0 && N < 2*m_step+1 ) { *real = m_FbRe[N]; *imag = m_FbIm[N]; }
};

void  EPG::GetZaState ( double* real, double* imag, const int &num ) const {
	if ( num >= 0 && num < m_step ) { *real = m_ZaRe[num]; *imag = m_ZaIm[num]; }
};

void  EPG::GetZbState ( double* real, double* imag, const int &num ) const {
	if ( num >= 0 && num < m_step ) { *real = m_ZbRe[num]; *imag = m_ZbIm[num]; }
};

/* return real and imaginary state values */
double EPG::GetReFa ( const int &num ) const { return m_FaRe[m_step-num-1]; }
double EPG::GetImFa ( const int &num ) const { return m_FaIm[m_step-num-1]; }
double EPG::GetReFb ( const int &num ) const { return m_FbRe[m_step-num];   } // the FB states shifted by 1, that means
double EPG::GetImFb ( const int &num ) const { return m_FbIm[m_step-num];   } // Fa_0 = Fa[m_step-1] ; Fb_0 = Fb[m_step]
double EPG::GetReZa ( const int &num ) const { return m_ZaRe[num]; }
double EPG::GetImZa ( const int &num ) const { return m_ZaIm[num]; }
double EPG::GetReZb ( const int &num ) const { return m_ZbRe[num]; }
double EPG::GetImZb ( const int &num ) const { return m_ZbIm[num]; }

/****************************************************************/
/* return magnitude state values */
double EPG::GetMagFa ( const int &num ) const {
		int N = m_step-num-1;
		if ( N >= 0 && N < 2*m_step+1 )	return sqrt (pow(m_FaRe[N],2) + pow(m_FaIm[N],2) ) ;
		else				return 0.0;
}
double EPG::GetMagZa ( const int &num ) const {
		if ( num >= 0 && num < m_step )	return sqrt (pow(m_ZaRe[num],2) + pow(m_ZaIm[num],2) ) ;
		else				return 0.0;
}
double EPG::GetMagFb ( const int &num ) const {
		int N = m_step-num;
		if ( N >= 0 && N < 2*m_step+1 )	return sqrt (pow(m_FbRe[N],2) + pow(m_FbIm[N],2) ) ;
		else				return 0.0;
}
double EPG::GetMagZb ( const int &num ) const {
		if ( num >= 0 && num < m_step )	return sqrt (pow(m_ZbRe[num],2) + pow(m_ZbIm[num],2) ) ;
		else				return 0.0;
}


double EPG::GetNextMagFa  ( const double &fa, const double &phi,const int &num  )  {
	m_step++;
  	Rotation(fa,phi,num,num);  // rotate single configuration
	double mag = GetMagFa (num);
  	Rotation(-fa,phi,num,num); // rotate it back
	m_step--;
	return mag;
};

/****************************************************************/
void EPG::Rotation ( const double &fa, const double &phi, const int &first, int last ) {

    int N = m_step-1;
    last = (last<0?N:last);

    if (first<0 || first>last || last>=m_step) { cout << "EPG::rotation error. States out of bounds.\n"; }

    double farad  =  fa*EPG_PI/180.0;
    double phirad = phi*EPG_PI/180.0;

    /*	These equations follow from matrix product a = A* b with complex magnetization vectors a,b after
	and before the pulse, respectively: b = [F(k) F*(-k) Z(k)]_{before},  a = [F(k) F*(-k) Z(k)]_{after}
	A: transition matrix (see EPG literature, e.g. Scheffler, Concepts Magn Reson 1999;11:291–304).  
    */
    for (int i=first;i<=last; ++i) {
	double a,b,c,d,e,f,g,h,hb,gb,ec,fc;
  
	a = cos(farad/2); a *= a; /* set coefficients of transition matrix */
	b = sin(farad/2); b *= b; c = sin(farad); d = cos(farad); e = sin(phirad);
	f = cos(phirad); g = sin(2.0*phirad); h = cos(2.0*phirad); hb = h*b; gb = g*b; ec = e*c; fc = f*c;    
      
	m_FaRe[N+i] = a*m_FbRe[N+i]  + hb*m_FbRe[N-i] + gb*m_FbIm[N-i] + ec*m_ZbRe[i] + fc*m_ZbIm[i];
	m_FaIm[N+i] = a*m_FbIm[N+i]  - hb*m_FbIm[N-i] + gb*m_FbRe[N-i] - fc*m_ZbRe[i] + ec*m_ZbIm[i];
	m_FaRe[N-i] = hb*m_FbRe[N+i] + gb*m_FbIm[N+i] + a*m_FbRe[N-i]  + ec*m_ZbRe[i] - fc*m_ZbIm[i];
	m_FaIm[N-i] = gb*m_FbRe[N+i] - hb*m_FbIm[N+i] + a*m_FbIm[N-i]  - fc*m_ZbRe[i] - ec*m_ZbIm[i];
	m_ZaRe[i] = (-ec*m_FbRe[N+i] + fc*m_FbIm[N+i] - ec*m_FbRe[N-i] + fc*m_FbIm[N-i] + 2.0*d*m_ZbRe[i])/2.0;
	m_ZaIm[i] = (-fc*m_FbRe[N+i] - ec*m_FbIm[N+i] + fc*m_FbRe[N-i] + ec*m_FbIm[N-i] + 2.0*d*m_ZbIm[i])/2.0;
    }
};


/****************************************************************/
void EPG::Step ( const double &fa, const double &phi, const bool &RFSpoil) {

	m_step++; //increase EPG step counter

	if (RFSpoil) {
		m_phase += m_step*phi;
		m_phase = fmod(m_phase,360.0);
	}
	else { m_phase = phi; } // remember phase in case needed for phase locking

	// extending state space, which shifts the Fb states: Fb[N] -> Fb[N+1]  (N = m_step-1)
	for (int i=0; i<2; i++) {
	   m_FaRe.push_back(0.0);
	   m_FaIm.push_back(0.0);
	   m_FbRe.push_back(0.0);
	   m_FbIm.push_back(0.0);
	}
	
	m_ZaRe.push_back(0.0);
	m_ZaIm.push_back(0.0);
	m_ZbRe.push_back(0.0);
	m_ZbIm.push_back(0.0);
	
	Rotation(fa,m_phase);

	// shifting and relaxation of states 
	// note that the center state of  m_FbRe is shifted one up with respect to m_FaRe
	int N = m_step-1;
	for (int i=0;i<m_step;++i) {
  		m_FbRe[N+i] = m_FaRe[N+i]*m_E2;	// evolution of dephasing transversal states (real part)
  		m_FbIm[N+i] = m_FaIm[N+i]*m_E2;	// evolution of dephasing transversal states (imaginary part)
  		m_FbRe[N-i] = m_FaRe[N-i]*m_E2;	// evolution of rephasing transversal states (real part)
  		m_FbIm[N-i] = m_FaIm[N-i]*m_E2;	// evolution of rephasing transversal states (imaginary part)
  		m_ZbRe[i] = m_ZaRe[i]*m_E1;	    // evolution of longitudinal states (real part)
  		m_ZbIm[i] = m_ZaIm[i]*m_E1;	    // evolution of longitudinal states (imaginary part)
	}
	m_ZbRe[0] += m_M0*(1.0-m_E1); 	    // recovery of longitudinal ground state 

};

/****************************************************************/
void EPG::Steps ( const double &fa, const double &phi, const int &steps, const bool &RFSpoil) {
	for(int i=0;i<steps; ++i) { Step(fa,phi,RFSpoil); }
};

void EPG::Steps ( const double* fa, const double &phi, const int &steps, const bool &RFSpoil) {
	for(int i=0;i<steps; ++i) { Step(fa[i],phi,RFSpoil); }
};

void EPG::Steps ( const double* fa, const double* phi, const int &steps) {
	for(int i=0;i<steps; ++i) { Step(fa[i],phi[i]); }
};

vector<double> EPG::GetMagTrain ( const std::vector<double> &fa, const std::vector<double> &phi) {

    std::vector<double> sig(fa.size());

	for(int i=0;i<fa.size(); ++i) {
        Step(fa[i],phi[i]);
        sig[i] = GetMagFa();
    }
        return sig;
};

int EPG::StepsToSS ( const double &fa, const double &Qphi, const double &tol ) {
	double F,Fold=-1.0,phi=0.0;
	for (int i=0;i<EPG_MAXSIZE;++i) {
		phi += i*Qphi;
		phi = fmod(phi,360.0);
		Step(fa,phi);
		F = GetMagFa();
		if (fabs((F-Fold)/m_M0)<tol) /*steady state reached*/ { return i; }
		Fold = F;
	}
	if (m_verbose) cout << "Warning EPG::GetSteadyState did not converge. Maximum number of possible states exceeded.\n";
	return -1;	
};


/****************************************************************/
bool EPG::FindFlipAngleTrain(const int &length, double* fa, const double* Ftarget, const double &reduce, const int &num, const double &tol ){

	EPG backup = *this;  //create a copy of current state
	double l=0.0,u=90.0; // lower/upper bound for bisectioning

	for (int i=0; i<EPG_MAXITER; ++i) {

		bool success = true;
		fa[0]=0.5*(l+u);
		Step(fa[0],0.0);
		double F0 = GetMagFa(num); //signal after first pulse
		for (int j=1;j<length;++j) {
			fa[j] = FindFlipAngle(0.0,90.0,F0*Ftarget[j],0.0,num,tol);
  			if ( fabs(GetNextMagFa(fa[j],0,num)-F0*Ftarget[j])>EPG_TOL  ) success = false;
			Step(fa[j],0.0);
		}
		if (success) { l=fa[0]; } else { u=fa[0]; }
		if ( fabs(u-l)<EPG_FA_TOL && success ) {  //found valid FA train
			if (reduce>0.0 && reduce < 1.0) { //reduce FAs by factor originally in fa[0]
				*this = backup;
				fa[0] *= reduce;
				Step(fa[0],0.0);
				F0 = GetMagFa(num); //signal after first pulse
				for (int j=1;j<length;++j) {
					fa[j] = FindFlipAngle(0.0,90.0,F0*Ftarget[j],0.0,num,tol);
  					if ( fabs(GetNextMagFa(fa[j],0,num)-F0*Ftarget[j])>EPG_TOL  ) success = false;
					Step(fa[j],0.0);
				}
			}
			return success;
		}
		*this = backup;
	}

	if (m_verbose) cout << "Warning EPG::FindFlipAngleTrain did not converge. Maximum number of iterations exceeded.\n";
	return false;
};

/****************************************************************/
double EPG::FindFlipAngle(const double &ax, const double &bx, const double &Ftarget, const double &phi, const int &num, const double &tol) {

  double a,b,c;				/* Abscissae, descr. see above	*/
  double fa;				/* f(a)				*/
  double fb;				/* f(b)				*/
  double fc;				/* f(c)				*/

  a = ax;  b = bx;
  fa = Ftarget - GetNextMagFa (a,phi,num);
  fb = Ftarget - GetNextMagFa (b,phi,num);
  
  c = a;   fc = fa;

  for(;;)		/* Main iteration loop	*/
  {
    double prev_step = b-a;		/* Distance from the last but one*/
					/* to the last approximation	*/
    double tol_act;			/* Actual tolerance		*/
    double p;      			/* Interpolation step is calcu- */
    double q;      			/* lated in the form p/q; divi- */
  					/* sion operations is delayed   */
 					/* until the last moment	*/
    double new_step;      		/* Step at this iteration       */
   
    if( fabs(fc) < fabs(fb) )
    {                         		/* Swap data for b to be the 	*/
	a = b;  b = c;  c = a;          /* best approximation		*/
	fa=fb;  fb=fc;  fc=fa;
    }
    tol_act = 2*EPG_EPSILON*fabs(b) + tol/2;
    new_step = (c-b)/2;

    if( fabs(new_step) <= tol_act || fb == (double)0 )
      return b;				/* Acceptable approx. is found	*/

    			/* Decide if the interpolation can be tried	*/
    if( fabs(prev_step) >= tol_act	/* If prev_step was large enough*/
	&& fabs(fa) > fabs(fb) )	/* and was in true direction,	*/
    {					/* Interpolatiom may be tried	*/
	double t1,cb,t2;
	cb = c-b;
	if( a==c )			/* If we have only two distinct	*/
	{				/* points linear interpolation 	*/
	  t1 = fb/fa;			/* can only be applied		*/
	  p = cb*t1;
	  q = 1.0 - t1;
 	}
	else				/* Quadric inverse interpolation*/
	{
	  q = fa/fc;  t1 = fb/fc;  t2 = fb/fa;
	  p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
	  q = (q-1.0) * (t1-1.0) * (t2-1.0);
	}
	if( p>(double)0 )		/* p was calculated with the op-*/
	  q = -q;			/* posite sign; make p positive	*/
	else				/* and assign possible minus to	*/
	  p = -p;			/* q				*/

	if( p < (0.75*cb*q-fabs(tol_act*q)/2)	/* If b+p/q falls in [b,c]*/
	    && p < fabs(prev_step*q/2) )	/* and isn't too large	*/
	  new_step = p/q;			/* it is accepted	*/
					/* If p/q is too large then the	*/
					/* bissection procedure can 	*/
					/* reduce [b,c] range to more	*/
					/* extent			*/
    }

    if( fabs(new_step) < tol_act ) {	/* Adjust the step to be not less*/
      if( new_step > (double)0 )	/* than tolerance		*/
	 new_step = tol_act;
      else
	new_step = -tol_act;
    }

    a = b;  fa = fb;			/* Save the previous approx.	*/
    b += new_step;			/* Do step to a new approxim.	*/  
    fb = Ftarget - GetNextMagFa (b,phi,num);
    
    if( (fb > 0 && fc > 0) || (fb < 0 && fc < 0) )
    {                 			/* Adjust c for it to have a sign*/
      c = a;  fc = fa;                  /* opposite to that of b	*/
    }
  }

};

