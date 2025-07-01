
# distutils: language=c++

import numpy as np
cimport numpy as np

# ctypedef int MINT
ctypedef np.int64_t MINT
ctypedef np.int16_t FFINT
ctypedef FFINT[:] FLAT
ctypedef FFINT[:,:] TABLE
ctypedef MINT[:,:] INDEXTABLE
ctypedef MINT[:] INDEX
ctypedef MINT[::1] INDEXFLAT
ctypedef FFINT[::1] FLATCONTIG
ctypedef FFINT[:,::1] TABLECONTIG

