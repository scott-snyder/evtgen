# Code review

HepForge uses Phabricator, which has a slightly different philosophy than
GitHub/GitLab with regard to code review.
It envisages that you will still commit to a feature branch in your local
repository but rather than pushing that branch to the central repository
and opening a pull/merge request where the review occurs, you instead
upload a diff to the Differential tool and conduct the review within that
tool.
Only once the review is completed there you (or someone with the
appropriate permissions) can "land" your changes onto the target branch
(usually master) in the central repository.
By default this will squash all the commits in your local feature branch
and "land" them as a single commit, but this behaviour can be modified (see
below).
In order to make this process easier, the
[Arcanist](https://secure.phabricator.com/book/phabricator/article/arcanist/)
command line tool is provided.
However, while documentation for this tool is improving it is still rather
sparse or entirely missing for some features, so I try here to provide some
guidance based on my own experimentation.

## Installing Arcanist

The only prerequisite that I know of is the PHP command line interface, for
which there is a minimum version requirement, which is currently 5.5.0

Otherwise you simply do:
```sh
git clone https://github.com/phacility/arcanist.git
```
then make sure that `arcanist/bin` is added to your `PATH` or set an alias
using the full path to the `arc` executable.

## Procedure for reviewees

For the most straightforward workflow you can do the following:
* Firstly, if you have not already done so, create a new task in Maniphest
  to describe what you're trying to do (fix a bug, implement a new feature,
  etc.) by following the instructions in
  [Working with Maniphest Tasks](Tasks.md)
* Create a new local branch based on the branch to which you want to eventually
  apply the changes, usually `master`
* Make your changes and commit as normal but **do not push the branch to
  HepForge**!
* Once you are ready to create a Differential revision, you can run the
  following commands from within your local git repository:
  * `arc which <target branch name>`

    This should tell you:
    * which remote repository it is going to interact with
    * which commits will be sent to Differential
    * whether a new revision will be created or if an existing one will be
      updated

    If these are all as expected then you can proceed with standard options
    for the rest of the commands.
    Otherwise, you will have to provide some extra options to the commands to
    override the behaviour - run `arc help command_name` to get further info
    (or see below for some of them).

  * `arc cover`

    This can tell you possible reviewers for the changes you've made.

  * `arc diff <target branch name>`

    This will create or update a Differential revision from your local
    changes.
    Probably you can use this without specifying options but you might want
    to use some of the following:
    * `--create` forces creation of a new revision
    * `--update <revision ID>` forces update of a specific revision
    * `--draft` create or update the revision in draft mode, which allows you
      to look over your changes in the web interface before requesting review
    * `--reviewers <usernames>` specify the reviewers
    * `--cc <usernames>` specify other users to be cc'ed

    You will be prompted to edit an automatically generated summary
    message.
    One useful edit would be to mention the associated Maniphest task so
    that the task and revision are automatically linked.

* Once the Differential revision is created, take a look at it in your web
  browser and, if necessary, modify its settings so that it is:
  * Visible To: Public
  * Editable By: EvtGen
  * Tags: EvtGen

* Once you are happy with the revision, you can select "Request Review" from
  the drop-down menu at the very bottom of the page, optionally adding a
  comment to the reviewers, and pressing "Publish revision"

* The review then proceeds via the web interface, with comments and replies
  and so on (see also the section below on 'Procedure for reviewers').
  If changes are requested then you can edit your code accordingly, commit,
  and then call `arc diff` again (probably also best to check that things
  will proceed as you expect by running `arc which` first) to update the
  revision.

* Once the revision is accepted by the reviewers it can be landed.
  If you have the appropriate permissions then you can do so, following the
  instructions below.
  Otherwise, one of the project admins will do so for you.

## Landing an accepted revision

If you already have the latest state of the revision in your local repository
(usually the case if you are the author of the revision) you can proceed
immediately to run `arc land`.
Otherwise, you first need to obtain the latest state by running `arc patch`.
Both commands are described below:

* `arc patch D12345`

  This will create a new branch in your local repository that applies the
  changes in the revision `D12345`, where you should substitute the number of
  the desired revision.
  Prior to running the command make sure you have the target branch of the
  revision (usually master) checked out and up to date.

* `arc land [ref]`

  This will publish the changes in the accepted revision onto a branch in the
  central repository.
  Sometimes this can be run with no options, although using `--preview` first
  is **always** a good idea.
  Sometimes it will be necessary to specify the tip commit or to use the
  `--revision` option, decribed below, in order to specify the commits to be
  landed.
  Sometimes it will be necessary to specify the remote and branch onto which
  to land, using the `--onto` and `--onto-remote` options described below.
  Here is a brief summary of the most useful options (run `arc help land` for
  more information):
  * `--preview` only prints the commits that will be published but does not
    make any local or remote modifications
  * `--hold` prints the commits that will be published and prepares the local
    modifications to be pushed but does not actually push anything
  * `--onto <branch name>` specifies the branch on the remote repository onto
    which the changes should be published (usually `master`)
  * `--onto-remote <remote name>` specifies the remote repository onto which
    the changes should be published (usually `origin`)
  * `--revision <revision-identifier>` Land a specific revision, rather than
    determining revisions automatically from the commits that are landing
  * `--strategy <squash/merge>` defines whether to squash local commits when
    publishing (the default behaviour) or whether to retain all local commits
    and create a merge commit

  When editing the commit message, bear in mind if you want to resolve and
  close the associated Task, in which case you can add `closes T123`.

## Procedure for reviewers

The review should consist of the obvious part, reading the code changes, and
also running some automated checks.
Depending on the outcome of each of these parts, you will then need to post
appropriate comments on the Differential revision and then mark the revision
either as being accepted or as requiring changes.
It may be necessary to iterate several times before accepting.

To run the automatic checks you should:
* Get a copy of the code in a new branch in your local repository using `arc
  patch`, described above
* Push that new branch into the CERN GitLab mirror repository
* The CI should then run - check the outcome of the pipeline
* Report any problems in the CI when commenting on the revision, along with any
  comments you have on the code itself

Other general points to check:

* The changes should be described in History.md, preferrably providing a link
  to the Maniphest task and/or Differential revision
* Applying clang-format (version 10.0.0) should result in no changes (should be
  automatically checked in the CERN GitLab CI - otherwise check manually)
* Names of classes, functions, variables, etc. should follow the existing conventions
* Any new files should have the EvtGen copyright notice at the top

