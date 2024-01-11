## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

devtools::check_rhub() generated 2 errors, 0 warnings, and 3 notes. The errors 
are apparently due to the test server not having the program JAGS installed. 
JAGS (https://mcmc-jags.sourceforge.net) is required for our package and this is 
indicated in our DESCRIPTION file.

We also ran spelling::spell_check_package(). This found a number of words, which
were all variable or parameter names, function names, package names, 
common abbreviations, technical terms, proper names, fragments of URLs, or 
journal names.

devtools::revdep() brought up nothing.
