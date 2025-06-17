# Using a .pxd file gives us a separate namespace for
# the C++ declarations. Using a .pxd file also allows
# us to reuse the declaration in multiple .pyx modules.
from libcpp cimport bool
from libcpp.vector cimport vector

cdef extern from "epg.h":
   
    cppclass EPG:
        EPG(double, double, double, double) except +
        EPG(EPG&) except +
        operator=(EPG) except +
        void SetParameters(double, double, double, double)
        void DeleteStates()

        int  GetStep()
        double GetM0()
        double GetT1()
        double GetT2()
        double GetTR()

        bool GetVerbose()
        void SetVerbose(bool)
        void Equilibrium()
        void NullTransverse()
        void SetStep(int)

        double GetMagFa(int)
        double GetMagFb(int)
        double GetMagZa(int)
        double GetMagZb(int)

        double GetReFa(int)
        double GetImFa(int)
        double GetReFb(int)
        double GetImFb(int)

        double GetReZa(int)
        double GetImZa(int)
        double GetReZb(int)
        double GetImZb(int)

        double GetNextMagFA(double, double, int)

        double GetPhase()

        void Step(double, double, bool)
        void Steps(double, double, int, bool)

        vector[double] GetMagTrain(vector[double] fa, vector[double] phi )

        int StepsToSS(double, double, double) ;

        bool   FindFlipAngleTrain(int, double* fa, double* Ftarget, double, int, double)

        double FindFlipAngle(double, double, double, double, int, double)
