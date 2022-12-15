# Making a release

When the master branch is ready to release there are several checks to perform
and, depending on the outcome of those checks, some resulting actions, and
finally there are further actions required to make the tag and tar file.


## Pre-release checks

The following items need to checked and, if found to need rectifying, the
appropriate changes should be made and committed to the master branch before
tagging it:

* Make sure that the History.md file contains appropriate entries for all
  changes since the previous tag
* Make sure that the History.md file contains the lines to indicate the new tag
* Make sure that the version number in the top-level CMakeLists.txt file
  matches that for the new tag
* Make sure that the version number in the setupEvtGen.sh script matches that
  for the new tag
* Double check that the version numbers in History.md, the top-level
  CMakeLists.txt, and setupEvtGen.sh all match


## Tagging the repository

To create the new tag, issue the command:
```sh
git tag -a RXX-YY-ZZ
```
where XX, YY, and ZZ should be replaced by the appropriate major, minor, and
patch version numbers, e.g. `R02-01-01`.
You will also need to explicitly push the tag:
```sh
git push RXX-YY-ZZ
```


## Create and upload the tar file

An easy way to create the tar file is to use the `git archive` command from
within your working copy, e.g.
```sh
git archive --prefix EvtGen/RXX-YY-ZZ/ -o EvtGen-XX.YY.ZZ.tar.gz RXX-YY-ZZ
```

The resulting file then has to be placed in the directory:
`/hepforge/projects/evtgen/downloads`
e.g. using `rsync`:
```sh
rsync -va EvtGen-XX.YY.ZZ.tar.gz login.hepforge.org:/hepforge/projects/evtgen/downloads/
```
NB the ssh port on login.hepforge.org is 222, so you should have an appropriate entry in your `~/.ssh/config` file, e.g.
```
Host login.hepforge.org
  User <your_username>
  Port 222
  PubkeyAuthentication yes
  ForwardX11 no
```

Make sure that the group ownership is `evtgen` and the permissions are `-rw-rw-r--.`.
To do this, ssh into `login.hepforge.org` and issue the commands:
```sh
chgrp evtgen /hepforge/projects/evtgen/downloads/EvtGen-XX.YY.ZZ.tar.gz
chmod 664 /hepforge/projects/evtgen/downloads/EvtGen-XX.YY.ZZ.tar.gz
```

