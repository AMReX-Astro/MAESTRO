# MAESTRO

> **_NOTE:_**  MAESTRO is no longer being actively developed.  Users should switch
> to MAESTROeX to take advantage of the latest capabilities: https://github.com/AMReX-Astro/MAESTROeX

*a low Mach number stellar hydrodynamics code*

`Maestro` solves the equations of low Mach number hydrodynamics for
stratified atmospheres/stars with a general equation of state.  It
includes reactions and thermal diffusion and can be used on anything
from a single core to 100,000s of processor cores with MPI + OpenMP.

A description of the algorithm and links to the algorithm papers can
be found here:

http://amrex-astro.github.io/MAESTRO/


## Getting Started

To use `Maestro` you need a copy of the ` FBoxLib` library, available
on github at:

https://github.com/AMReX-Codes/FBoxLib.git

To use anything other than the simple microphysics, you need the
StarKiller microphysics package, available on github at:

https://github.com/starkiller-astro/Microphysics.git

There are a few environment variables that need to be set.  A *Getting
Started* guide is provided in the Maestro User's Guide which will walk
you through this.  To build the User's Guide, cd into `Docs/`, and
type `make`, or download the PDF here:

http://bender.astro.sunysb.edu/Maestro/staging/MAESTRO/Docs/MaestroUsersGuide.pdf


## Call Graph

A doxygen-generated call graph is available here:

http://bender.astro.sunysb.edu/Maestro/staging/MAESTRO/html/


## Development Model

Development generally follows the following ideas:

  * New features are committed to the `development` branch.

    Nightly regression testing is used to ensure that no answers
    change (or if they do, that the changes were expected).

    If a change is critical, we can cherry-pick the commit from
    `development` to `master`.

  * Contributions are welcomed from anyone.  *Any contributions that
    have the potential to change answers should be done via pull
    requests.*   A pull request should be generated from your fork of
    Maestro and target the `development` branch.  (If you mistakenly
    target `master`, we can change it for you.)

    Please add a line to `CHANGES` summarizing your change if it
    is a bug fix or new feature.  Reference the PR or issue as
    appropriate

    If there are a number of small commits making up the PR, we may
    wish to squash commits upon merge to have a clean history.
    *Please ensure that your PR title and first post are descriptive,
    since these will be used for a squashed commit message.*

  * On the first workday of each month, we perform a merge of
    `development` into `master`, in coordination with `AMReX` and
    `FBoxLib`, `Maestro`, and `Microphysics`.  For this merge to take
    place, we need to be passing the regression tests.

    To accommodate this need, we close the merge window into
    `development` a few days before the merge day.  While the merge
    window is closed, only bug fixes should be pushed into
    `development`.  Once the merge from `development` -> `master` is
    done, the merge window reopens.


## Core Developers

People who make a number of substantive contributions will be named
"core developers" of Maestro.  The criteria for becoming a core
developer are flexible, but generally involve one of the following:

  * 10 non-merge commits to `MAESTRO/Source/` or `MAESTRO/Docs/` or one
    of the problems that is not your own science problem *or*

  * addition of a new algorithm / module  *or*

  * substantial input into the code design process or testing

Core developers will be recognized in the following ways:

  * invited to the group's slack team

  * listed in the User's Guide and website as a core developer

  * invited to co-author general code papers / proceedings describing
    Maestro, its performance, etc.  (Note: science papers will always
    be left to the science leads to determine authorship).

If a core developer is inactive for 3 years, we may reassess their
status as a core developer.


## Mailing list:

You can subscribe to the maestro-help mailing list at google groups at:
https://groups.google.com/forum/#!forum/maestro-help




