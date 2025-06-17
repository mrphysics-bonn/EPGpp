# distutils: language = c++
# distutils: sources = epg.cpp

from epg cimport EPG
cimport numpy as np
import numpy as np
from libcpp.vector cimport vector
from libc.stdlib cimport free
from cpython cimport PyObject, Py_INCREF
from cython.operator cimport dereference as deref

# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()

#FIXME: Get constants from epg.h directly
cdef EPG_TOL = 0.0001

cdef extern from "<utility>" namespace "std" nogil:
    vector[double] move(vector[double])

# define ArrayWrapper as holding in a vector
cdef class ArrayWrapper:
    cdef vector[double] vec
    cdef Py_ssize_t shape[1]
    cdef Py_ssize_t strides[1]

    # constructor and destructor are fairly unimportant now since
    # vec will be destroyed automatically.

    cdef set_data(self, vector[double]& data):
       #self.vec = move(data)
       self.vec.swap(data)
       # @ead suggests `self.vec.swap(data)` instead
       # to avoid having to wrap move

    # now implement the buffer protocol for the class
    # which makes it generally useful to anything that expects an array
    def __getbuffer__(self, Py_buffer *buffer, int flags):
        # relevant documentation http://cython.readthedocs.io/en/latest/src/userguide/buffer.html#a-matrix-class
        cdef Py_ssize_t itemsize = sizeof(self.vec[0])

        self.shape[0] = self.vec.size()
        self.strides[0] = sizeof(double)
        buffer.buf = <char *>&(self.vec[0])
        buffer.format = 'd'
        buffer.internal = NULL
        buffer.itemsize = itemsize
        buffer.len = self.vec.size() * itemsize   # product(shape) * itemsize
        buffer.ndim = 1
        buffer.obj = self
        buffer.readonly = 0
        buffer.shape = self.shape
        buffer.strides = self.strides
        buffer.suboffsets = NULL

cdef class PyEPG:
    """
    Cython wrapper class for C++ class EPG

    ...

    Attributes
    ----------
    M0 : double
        Equilibrium magnetization [a.u.]
    T1 : double
        Longitudinal relaxation time [ms]
    T2 : double
        Transverse relaxation time [ms]
    TR : double
        repetition time [ms]

    Methods
    -------
    SetParameters(M0, T1, T2,TR)
        Sets simulation parameters
    DeleteStates()
        Clears all state vectors (Fa, Fb, Za, Zb)
    GetStep()
        Return current step of this EPG
    GetM0()
        Return equilibrium magnetization
    GetT1()
        Return longitudinal relaxation
    GetT2()
        Return transverse relaxation
    Equilibrium()
        Reset EPG to equilibrium
    NullTransverse()
        Set all transverse state to zero (aKa spoiling)
    SetStep(val):
        Sets step counter of this EPG
    GetMagFa(int num=0):
        Returns transverse magnitude of state num after last RF pulse
    GetMagFb(int num=0):
        Returns transverse magnitude of state num before next RF pulse
    GetMagZa(int num=0):
        Returns longitudinal magnitude of state num after last RF pulse
    GetMagZb(int num=0):
        Returns longitudinal magnitude of state num before next RF pulse
    Step(fa, phi, RFSpoil=true):
        Single EPG step forward
    Steps(fa, phi, TR, steps):
        Multiple EPG steps forward, constant flipangle and phase
    GetMagTrain(fa, phi, TR):
        Multiple EPG steps forward, variable flipangle and phase
    StepsToSS(fa, Qphi, TR, tol=EPG_TOL):
        Stepping until steady state is reached
    FindFlipAngleTrain(np.ndarray[double, ndim=1, mode="c"] Ftarget, TR, creduce=0.0, num=0, tol=EPG_TOL):
        Find fli pangles for a given target signal shape
    """

    cdef EPG* _thisptr

    def __cinit__(self, other=1.0, T1=1000.0, T2=100.0, TR=10.0):
        # Initialize the "this pointer" to NULL so __dealloc__
        # knows if there is something to deallocate. Do not
        # call new EPG() here.
        # self._thisptr = NULL
        cdef PyEPG ostr
        if other and type(other) is PyEPG:
            ostr = <PyEPG> other
            self._thisptr = new EPG(deref(ostr._thisptr))
        else:
            self._thisptr = new EPG(other, T1, T2, TR)

    def __dealloc__(self):
        # Only call del if the C++ object is alive,
        # or we will get a segfault.
        if self._thisptr != NULL:
            del self._thisptr

    cdef int _check_alive(self) except -1:
        # Beacuse of the context manager protocol, the C++ object
        # might die before PyEPG self is reclaimed.
        # We therefore need a small utility to check for the
        # availability of self._thisptr
        if self._thisptr == NULL:
            raise RuntimeError("Wrapped C++ object is deleted")
        else:
            return 0


 # EPG state creation and initialization
    def SetParameters(PyEPG self, M0, T1, T2, TR):
        """
        Sets simulation parameters


        Parameters
        ----------
        M0 : double
            Equilibrium magnetization [a.u.]
        T1 : double
            Longitudinal relaxation time [ms]
        T2 : double
            Transverse relaxation time [ms]
        TR : double
            repetition time [ms]
        """

        self._check_alive()
        self._thisptr.SetParameters(M0, T1, T2, TR)

    def DeleteStates(PyEPG self):
        """
        Clears all state vectors (Fa, Fb, Za, Zb), leaving the vectors with a size of 0.
        """

        self._check_alive()
        return self._thisptr.DeleteStates()

    def GetStep(PyEPG self):
        """
        Return current step of this EPG
        """

        self._check_alive()
        return self._thisptr.GetStep()

    def GetM0(PyEPG self):
        """
        Return equilibrium magnetization
        """

        self._check_alive()
        return self._thisptr.GetM0()

    def GetT1(PyEPG self):
        """
        Return longitudinal relaxation
        """

        self._check_alive()
        return self._thisptr.GetT1()

    def GetT2(PyEPG self):
        """
        Return transverse relaxation
        """

        self._check_alive()
        return self._thisptr.GetT2()

    def GetVerbose(PyEPG self):
        self._check_alive()
        return self._thisptr.GetVerbose()

    def SetVerbose(PyEPG self, val):
        self._check_alive()
        self._thisptr.SetVerbose(val)

    def Equilibrium(PyEPG self):
        """
        Reset EPG to equilibrium
        """

        self._check_alive()
        self._thisptr.Equilibrium()

    def NullTransverse(PyEPG self):
        """
        Set all transverse state to zero (aKa spoiling)
        """

        self._check_alive()
        self._thisptr.NullTransverse()

    def SetStep(PyEPG self, val):
        """
        Sets step counter of this EPG


        Parameters
        ----------
        val : int
            step counter
        """

        self._check_alive()
        self._thisptr.SetStep(val)

    def GetMagFa(PyEPG self, int num=0):
        self._check_alive()
        return self._thisptr.GetMagFa(num)

    def GetMagFb(PyEPG self, int num=0):
        self._check_alive()
        return self._thisptr.GetMagFb(num)

    def GetMagZa(PyEPG self, int num=0):
        self._check_alive()
        return self._thisptr.GetMagZa(num)

    def GetMagZb(PyEPG self, int num=0):
        self._check_alive()
        return self._thisptr.GetMagZb(num)

    def GetReFa(PyEPG self, int num=0):
        self._check_alive()
        return self._thisptr.GetReFa(num)

    def GetImFa(PyEPG self, int num=0):
        self._check_alive()
        return self._thisptr.GetImFa(num)

    def GetReFb(PyEPG self, int num=0):
        self._check_alive()
        return self._thisptr.GetReFb(num)

    def GetImFb(PyEPG self, int num=0):
        self._check_alive()
        return self._thisptr.GetImFb(num)

    def GetReZa(PyEPG self, int num=0):
        self._check_alive()
        return self._thisptr.GetReZa(num)

    def GetImZa(PyEPG self, int num=0):
        self._check_alive()
        return self._thisptr.GetImZa(num)

    def GetReZb(PyEPG self, int num=0):
        self._check_alive()
        return self._thisptr.GetReZb(num)

    def GetImZb(PyEPG self, int num=0):
        self._check_alive()
        return self._thisptr.GetImZb(num)

    def GetPhase(PyEPG self):
        self._check_alive()
        return self._thisptr.GetPhase()

    def Step(PyEPG self, fa, phi, RFSpoil = False):
        """
        Single EPG step forward

        Parameters
        ----------
        fa : double
            RF flipangle [deg]
        phi : double
            RF phase
        RFSpoil: bool 
            RFSpoiling yes/no? default: False
        
        Remarks
        -------
        1 Step = RF-pulse + time evolution in TR

        """

        self._check_alive()
        self._thisptr.Step(fa, phi, RFSpoil)

    def Steps(PyEPG self, fa, phi, steps, RFSpoil = False):
        """
        Multiple EPG steps forward, constant flipangle and phase

        Parameters
        ----------
        fa : double
            RF flipangle [deg]
        phi : double
            RF phase
        steps : int
            number of steps
        RFSpoil: bool 
            RFSpoiling yes/no? default: False

        Remarks
        -------
        n Steps = n * (RF-pulse + TR intervall)

        """

        self._check_alive()
        self._thisptr.Steps(fa, phi, steps, RFSpoil)

    def GetMagTrain(PyEPG self, fa, phi):
        """
        Multiple EPG steps forward, variable flipangle and phase

        Parameters
        ----------
        fa : 1d-array (double)
            RF flipangle [deg]
        phi : 1d-array (double)
            RF phase
        steps : int
            number of steps

        Returns
        -------
        signal : 1d-array (double)
            Transverse magnitude of state zero

        """
        self._check_alive()
        cdef vector[double] array = self._thisptr.GetMagTrain(<vector[double]&> fa, <vector[double]&> phi)

        w = ArrayWrapper()
        w.set_data(array) # "array" itself is invalid from here on
        ndarray = np.asarray(w)

        return ndarray

    def StepsToSS(PyEPG self, fa, Qphi, tol=EPG_TOL):
        """
        Stepping until steady state is reached


        Parameters
        ----------
        fa : double
            RF flipangle [deg]
        Qphi : double
            Quadratic RF phase increment
        tol : double
           tolerance for termination

        Returns
        -------
        steps : 1d-array (double)
            Steps needed to reach steady state

        """
        self._check_alive()
        return self._thisptr.StepsToSS(fa, Qphi, tol)

    def FindFlipAngleTrain(PyEPG self, np.ndarray[double, ndim=1, mode="c"] Ftarget, creduce=0.0, num=0, tol=EPG_TOL):
        """
        Find flip angles for a given target signal shape

        Parameters
        ----------
        Ftarget : 1d-array (double)
            Target signal shape
        creduce : double
            reduction factor for flip angle train scaling (0.0 < creduce < 1.0)
        num : int
            EPG state order (defaults to zero)
        tol : double
           tolerance for termination

        """
        # Good example how to deal with C-pointers and Numpy arrays
        self._check_alive()

        length = Ftarget.size
        fa = np.zeros(length, dtype=np.double)

        cdef np.ndarray[double, ndim=1, mode="c"] np_buff = np.ascontiguousarray(fa, dtype=np.double)

        self._thisptr.FindFlipAngleTrain(length, &np_buff[0], &Ftarget[0], creduce, num, tol)
        return np_buff

    # The context manager protocol allows us to precisely
    # control the liftetime of the wrapped C++ object. del
    # is called deterministically and independently of
    # the Python garbage collection.

    def __enter__(PyEPG self):
        self._check_alive()
        return self

    def __exit__(PyEPG self, exc_tp, exc_val, exc_tb):
        if self._thisptr != NULL:
            del self._thisptr
            self._thisptr = NULL # inform __dealloc__
        return False # propagate exceptions

    def FindFlipAngle(PyEPG self,  FAmin, FAmax, FMagTarget, phi=0, num=0, tol=EPG_TOL): 
        """
        Find flip angle within bounds to match a target magnitude signal in the next step

        Parameters
        ----------
        FAmin : double
            lower bound = minimum RF flipangle [deg]
        FAmin : double
            upper bound = maximum RF flipangle [deg]
        FMagTarget : double
            target magnitude signal
        phi : double
           RF pulse phase
        num : int
           EPG state order (defaults to zero)
        tol : double
           tolerance for termination

        Returns
        -------
        FA : double
            flip angle which generates target signal (might fail)

        """
        self._check_alive()
        return self._thisptr.FindFlipAngle(FAmin, FAmax, FMagTarget, phi, num, tol)