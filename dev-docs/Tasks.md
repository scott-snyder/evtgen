# Working with Maniphest Tasks

In the HepForge Phabricator instance, tasks are managed by the
[Maniphest](https://phab.hepforge.org/maniphest/) application.
The tasks associated with EvtGen are collected in the
[EvtGen Workboard](https://phab.hepforge.org/project/board/135/), where
they are organised into columns depending on the task category (bug
reports, new features or models, modernisation, documentation).

To create a new task you can click on the arrow in the top right corner of
the column in which the task belongs and select "Create task".
In the box that pops up you need to supply:
* A title, which should be a very short summary of the task
* Optionally you can assign a user to the task
* Optionally you can set the priority
* Supply a more detailed description
* Set the visibility to Public
* Set the editability to EvtGen
* Set the tags to EvtGen
* Optionally add subscribers, e.g. evtgen Admins

It is also now possible to quickly create a task (which can later be
fleshed-out using the web interface) using the
[Arcanist](https://secure.phabricator.com/book/phabricator/article/arcanist/)
tool (see the [code review](CodeReview.md) documentation for instructions
on installing Arcanist) by doing:
```sh
arc todo "Title of task" --project EvtGen
```
Adding the `--browse` option will automatically open the newly created task
in a web browser.

