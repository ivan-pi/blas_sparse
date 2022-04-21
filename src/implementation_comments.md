Here are a few points concerning the implementation of the Fortran interface:

* The F95 interface can be implemented as a thin wrapper of the C interface
* With some caveats, it could also be implemented as wrapper of the F77 interface, but not in a standard-conforming way

Some of the problematic areas include:
* assumed-shape arrays can be discontiguous, i.e. strided arrays. This raises the issue how to correctly determine stride, and prevent unwanted temporary copies
* the Fortran interface might diverge from the C implementation in error handling, e.g. we might want to add status flags which check conformity of the input arguments
* assuming we want to minimize the number of indirections, we need to either
  - call the C procedures in Fortran; this is not possible for all due to the stride issue
  - implement the procedures directly in C; this requires bind(C) on the Fortran side which has some caveats; we also may need to use "ISO_Fortran_binding.h"

There are more options of course
- implement the procedures in Fortran directly. This would be the most desirable in some ways
- implement the procedures in C, but directly upon the internal array structure; since this involves quite a lot of pointer trickery, it's not for the faint-hearted
- A kind of middle ground, would be to actually implement the algorithm in C++, but upon containers of Fortran