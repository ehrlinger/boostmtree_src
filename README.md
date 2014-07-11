boostmtree
========================================================================

# The current state of your session:
> git config --list

# If git doesn't know who you are:
> git config --global user.email "commerce@kogalur.com"
> git config --global user.name "kogalur"

# Cloning a repository:

> cd ~/boosmtree/cloud (or where ever you want it)

> git clone https://github.com/kogalur/boostmtree/

#You will also want to make a user friendly name for each cloned directory, since that is all local:

# User friendly name for remote:
> git remote add boostmtree https://github.com/kogalur/boostmtree/

# user name and password access
UserID:       kogalur
Password:     MW4euZE8nP%Ndx


If you've added new files:
> git add cran/install/R/new.file.R

Finally, when you're ready to post the build, you must first commit.
This tells the cloud what you intend to push, but does not actually
push it to the cloud.  The dot is just the associated message, but we
have been relying on the README file for more annotated information:
> git commit -am "."

Then, we push the changes up to the cloud:
> git push boostmtree

Then we want to tag the build so we remember the snapshot:
> git tag -a bld2014xxxx -m '.'

The dot is just the associated message, but we have been relying on the README
file for more annotated information.
> git tag -a bld2014xxxx -m "."

Or, if you need to move a tag or retag after minor corrections:
> git tag -f -a bld2014xxxx -m "."

Then you have to push the tag to the cloud:
> git push boostmtree bld2014xxxx

To view all tags
> git tag

To review the details of a tag
> git show bld2014xxxx



# Pulling from the trunk:
> git pull boostmree master


## NICETIES:

# Viewing remote aliases:
> git remote -v

# User friendly name for remote:
> git remote add boostmtree https://github.com/kogalur/boostmtree/

# Check that your alias is in effect:
> git remote -v

# Viewing changes:
> git log --stat --summary

gives a summary of changes of all commits.

> git log -p

shows the diffs at each commit.


------------------------------------------------------------------------
------------------------------------------------------------------------
TAG             INIT    DESCRIPTION
------------------------------------------------------------------------
------------------------------------------------------------------------
bld20140711     ubk     Revert to bld20140616
bld20140710     hi      Bad commit by Hemant.
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
