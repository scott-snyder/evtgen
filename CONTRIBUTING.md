# Guidelines for contributing bug reports or code to EvtGen

In order to submit bug reports or to submit new/modified code for review, you
will need an account on [HepForge](https://www.hepforge.org/register).

Once your account is created, you should contact the
[EvtGen developers](mailto:evtgen@projects.hepforge.org) to request that your
account be added to the EvtGen project.

## Submitting bug reports

In the HepForge Phabricator instance, bug reports are handled by the
[Maniphest](https://phab.hepforge.org/maniphest/) application.
Please follow the instructions in [Working with Maniphest Tasks](dev-docs/Tasks.md),
making sure to use the *Bug reports* column on the
[EvtGen Workboard](https://phab.hepforge.org/project/board/135/).

## C++ standard and coding conventions

New EvtGen code can assume the use of the C++ 17 standard.

We will be introducing certain static analysis checks (using `clang-tidy`).
These will be documented here as they are introduced.
Almost certainly they will include modernisation checks such as preferring use of:
* `nullptr`
* uniform initialisation
* override
* etc.

As such, any new code should try to conform to these already.

## Code formatting, naming conventions, etc.

EvtGen uses a `clang-format` style file for code formatting.
Please run `clang-format` version 16 on all files that you have modified/created.

EvtGen is licenced under the GNU GPL v3 (or later) and as such the appropriate
copyright header should appear at the top of each source file.

The naming conventions in the code have not been very strictly applied in the
past and so there is currently some variation throughout the source code.
We are attempting to gradually unify the naming of all classes, functions,
variables, etc. and update the guidelines here accordingly.
The current guidelines are:
* All class/struct/namespace/enumeration names should be UpperCamelCase and
  should start with the Evt prefix
* All enumeration states should also be UpperCamelCase
* All function and variable names should be lowerCamelCase
* All member variables should be prefixed with `m_`
* When including headers, only includes of the C++ standard library headers should use the `#include <...>` form, all others should use the `#include "..."` form

As specified in [T105](https://phab.hepforge.org/T105), the use of inclusive
language is preferred and encouraged in the EvtGen project.
New code is encouraged to adopt gender-neutral terminology as much as possible.
Existing gendered terminology will gradually be migrated to more neutral terms.

## Additional checks before submitting code for review

* Have you documented your changes in the History.md file, referencing the
  appropriate Maniphest task and/or Differential revision?

## Submitting code for review

Follow the documented procedures for [code review](dev-docs/CodeReview.md).
