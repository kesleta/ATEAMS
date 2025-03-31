
# distutils: language=c++
import cython
from libcpp.unordered_map cimport unordered_map as Map
cdef Map[int, int[]] P;
