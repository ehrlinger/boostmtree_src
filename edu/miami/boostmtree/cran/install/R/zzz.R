.onAttach <- function(libname, pkgname) {
    boostmtree.version <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                      fields="Version")
    packageStartupMessage(paste("\n",
                                pkgname,
                                boostmtree.version,
                                "\n",
                                "\n",
                                "Type boostmtree.news() to see new features, changes, and bug fixes.",
                                "\n",
                                "\n"))
}
