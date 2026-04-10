boostmtree.news <- function(...) {
    newsfile <- file.path(system.file(package="boostmtree"), "NEWS.md")
    file.show(newsfile)
}
