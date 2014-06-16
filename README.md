boostmtree
========================================================================


If you've added new files:

> git add cran/install/R/new.file.R

Finally, when you're ready to post the build, you first commit:

> git commit -am "."

This tells the cloud what you intend to push.  Then tagging, so we remember the snapshot:

> git tag -am "."  bld2014xxxx

Or, if you need to move a tag or retag

> git tag -f -am "." bld2014xxxx


Then, finally we push up to the cloud:

> git push boostmtree


# Viewing remote aliases:
git remote -v

# User friendly name for remote:
git remote add boostmtree https://github.com/kogalur/boostmtree/

# Viewing remote aliases after you have added the above:
git remote -v

# Pulling from the trunk:
git pull boostmree master

# Cloning a repository:

> cd ~/boosmtree/cloud (or where ever you want it)

> git clone https://github.com/kogalur/boostmtree/

#You will also want to make a user friendly name for each cloned directory, since that is all local:

# User friendly name for remote:
> git remote add boostmtree https://github.com/kogalur/boostmtree/

# user name and password access

UserID:           kogalur
Password:     MW4euZE8nP%Ndx

------------------------------------------------------------------------
------------------------------------------------------------------------
TAG             INIT    DESCRIPTION
------------------------------------------------------------------------
------------------------------------------------------------------------
bld20140616     hi      Compiled man pages with illustrative examples.
                        Upgraded boostmtree code.  Y-values are scaled.
                        Added forest base-learners (TBD predict mode).
                        Adaptive smoothing improved by using the fixed
                        lambda gamma update in place of BLUP.  Resolved
                        long standing issues with phi and rho by noting
                        phi is a nuisance parameter. Removed the M
                        burn-in option.  Added hidden variables ntree,
                        bootstrap and rho.  Added several new functions 
                        to utilities including hidden variable parsing.  
                        Predict code upgraded to accomodate scaling of Y.  
                        Replaced benchmark test code with improved examples. 

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
