boostmtree
========================================================================

# Create a directory for your project called "boostmtree"
# in your user directory

> mkdir ~/working
> cd ~/working

# Set up and initialize the necessary Git files in ~/working/.git
> git init

# Cloning the repository in the cloud to your local machine.
> git clone https://github.com/kogalur/boostmtree/

# Move to the source tree "home" directory
> cd boostmtree/edu/miami/boostmtree

# Adding files
> git add cran/install/R/new.file.R

# Commiting files
> git commit -am "."

# Tagging:
> git tag -a bld2014xxxx
> git tag -f -am "." bld2014xxxx

# Viewing remote aliases:
git remote -v

# User friendly name for remote:
git remote add boostmtree https://github.com/kogalur/boostmtree/

# Viewing remote aliases after you have added the above:
git remote -v

# Pushing to the cloud:
git push boostmtree

# Pulling from the trunk:
git pull boostmree master


------------------------------------------------------------------------
------------------------------------------------------------------------
TAG             INIT    DESCRIPTION
------------------------------------------------------------------------
------------------------------------------------------------------------
bld20140615     hi      Compiled man pages with illustrative examples.
                        Upgraded boostmtree code.  Y-values are scaled.
                        Adpative smoothing improved by using the fixed
                        lambda gamma update in place of BLUP.  Resolved
                        long standing issues with phi and rho by noting
                        phi is a nuisance parameter. Removed the M
                        burn-in option.  Added hidden variables ntree,
                        bootstrap and rho.  Added several new functions 
                        to utilities.  Predict code upgraded to accomodate 
                        scaling of Y.  Replaced benchmark test code with 
			improved examples. 

------------------------------------------------------------------------
TRUNK TAGS
------------------------------------------------------------------------
bld20140610     ubk     Initial check in of all code and build framework


------------------------------------------------------------------------
RELEASE TAGS
------------------------------------------------------------------------
release_1_0_0   ubk     Equivalent to tag bld2014xxxx


------------------------------------------------------------------------
BRANCH TAGS
------------------------------------------------------------------------

Blah blah blah blah
------------------------------------------------------------------------
bra2014xxxx     ubk         Blah blah blah blah
