# MAESTRO

*a low Mach number stellar hydrodynamics code*

`Maestro` solves the equations of low Mach number hydrodynamics for
stratified atmospheres/stars with a general equation of state.  It
includes reactions and thermal diffusion and can be used on anything
from a single core to 100,000s of processor cores with MPI + OpenMP.

A description of the algorithm and links to the algorithm papers can
be found here:

http://boxlib-codes.github.io/MAESTRO/


## Getting Started:

To use `Maestro` you need a copy of the `BoxLib` library, available
on github at:

https://github.com/BoxLib-Codes/BoxLib.git

A *Getting Started* guide is provided in the Maestro User's Guide.  To
build the User's Guide, cd into `Docs/`, and type `make`, or download
the PDF here:

http://bender.astro.sunysb.edu/Maestro/staging/MAESTRO/Docs/MaestroUsersGuide.pdf


## Call Tree:

doxygen-generated call trees are available here:

http://bender.astro.sunysb.edu/Maestro/staging/MAESTRO/html/


## Development Model:

New features are committed to the `development` branch.  Nightly
regression testing is used to ensure that no answers change (or if
they do, that the changes were expected).  No changes should ever
be pushed directly into `master`.

On the first workday of each month, we perform a merge of
`development` into `master`, in coordination with `BoxLib`, `Castro`,
and `Microphysics`.  For this merge to take place, we need to be
passing the regression tests.  To accommodate this need, we close the
merge window into `development` a few days before the merge day.
While the merge window is closed, only bug fixes should be pushed into
`development`.  Once the merge from `development` -> `master` is done,
the merge window reopens.


## Mailing list:

You can subscribe to the maestro-help mailing list at google groups at:
https://groups.google.com/forum/#!forum/maestro-help




