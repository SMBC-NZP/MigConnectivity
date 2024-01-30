# Prepare for CRAN ----

# Update dependencies in DESCRIPTION
# install.packages('attachment', repos = 'https://thinkr-open.r-universe.dev')
attachment::att_amend_desc()

# Run tests and examples
devtools::test()
devtools::run_examples()
# autotest::autotest_package(test = TRUE)

# Check package as CRAN
results <- rcmdcheck::rcmdcheck(args = c("--no-manual", "--as-cran"))
results
# Check content
# install.packages('checkhelper', repos = 'https://thinkr-open.r-universe.dev')
checkhelper::find_missing_tags()
# _Check that you left the house clean after the check, examples and tests
all_files_remaining <- checkhelper::check_clean_userspace()
print(all_files_remaining, n = nrow(all_files_remaining))
write.csv(all_files_remaining, "data-raw/dev/files_remaining_log.csv")

# Check spelling
# usethis::use_spell_check()
spelling::spell_check_package()

# Check URL are correct
# install.packages('urlchecker', repos = 'https://r-lib.r-universe.dev')
urlchecker::url_check(parallel = F)
urlchecker::url_update()

# check on other distributions
# _rhub
theres_the_rhub <- devtools::check_rhub()
print(theres_the_rhub)
rhub::check_on_windows(check_args = "--force-multiarch")
rhub::check_on_solaris()
# _win devel CRAN
devtools::check_win_devel()
# _macos CRAN
devtools::check_mac_release()

# Check reverse dependencies
# remotes::install_github("r-lib/revdepcheck")
install.packages('revdepcheck', repos = 'https://r-lib.r-universe.dev')
usethis::use_git_ignore("revdep/")
usethis::use_build_ignore("revdep/")

devtools::revdep()
library(revdepcheck)
# In another session
id <- rstudioapi::terminalExecute("Rscript -e 'revdepcheck::revdep_check(num_workers = 4)'")
rstudioapi::terminalKill(id)
# See outputs
revdep_details(revdep = "pkg")
revdep_summary()                 # table of results by package
revdep_report() # in revdep/
# Clean up when on CRAN
revdep_reset()

# Update NEWS
# Bump version manually and add list of changes

# Add comments for CRAN
usethis::use_cran_comments(open = rlang::is_interactive())

# Upgrade version number
usethis::use_version(which = c("patch", "minor", "major", "dev")[1])

tools::buildVignettes(dir = ".", tangle=TRUE)
dir.create("docs")
file.copy(dir("vignettes", full.names=TRUE), "docs", overwrite=TRUE)

# Verify you're ready for release, and release
devtools::release()
devtools::submit_cran()
