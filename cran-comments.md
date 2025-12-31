## Test environments
* local OS X install, R 4.3.2
* win-builder (devel and release)

## R CMD check results
0 errors | 0 warnings | 0 notes

## Reverse dependencies
This is a new package, so there are no reverse dependencies.

## Resubmission
This is a resubmission. In this version I have:
* Removed unnecessary whitespace from the Description field.
* Added references (Tibshirani 1996, Meinshausen 2010) to the Description field.
* Added \value tags to .Rd files for exported functions (gradLasso, simulate_data).
* Added a 'verbose' argument to gradLasso() and cv.gradLasso() to suppress console output by default.
* Ensured user settings (par) are reset using on.exit() in plot methods.
