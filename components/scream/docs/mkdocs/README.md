# SCREAM Documentation

This "meta-documentation" is for you, the reader interested in contributing to
our nascent documentation effort. Welcome!

We are trying out [MkDocs](https://mkdocs.org), a simple static website generator
for people who are tired of trying to get Sphinx to work properly. Our goal is
to use MkDocs to build two important pieces of documentation for SCREAM:

1. A User Guide, which explains to climate scientists how to build SCREAM and
   configure and run cases of interest.
2. A Developer Guide, which serves as the primary resource for developers who
   contribute to SCREAM.

This documentation lives on [SCREAM's GitHub project page](https://e3sm-project.github.io/scream/).
Take a look!

You can read all you want about MkDocs at the above link, but here are the
essentials:

* In the MkDocs universe, documentation is written using Markdown and a set
  of nifty configurable extensions (including Latex-style math!).
* Documentation files are organized into folders within this one:
  * The User Guide, which explains how to configure and run SCREAM, lives
    in `user/`.
  * The Developer Guide, which helps SCREAM developers find out what's where,
    lives in `dev/`.
  * Documentation common to these two guides is stored in `common/`.
  In principle, we can make as deep a hierarchy as we want (and people with
  advanced degrees love to go nuts with deep hierachies!), but let's start
  simple.
* MkDocs is configured using settings in `mkdocs.yml`, which lives in SCREAM's
  top-level source directory.
* Typically, you won't have to edit `mkdocs.yml`. The only thing you'll
  probably have to modify is the nav section, near the top of the file. This
  section contains three basic elements:
  * The [SCREAM documentation home page](https://e3sm-project.github.io/scream/) itself
  * A link to the User Guide and its sections
  * A link to the Developer Guide and its sections

## How Do I Use MkDocs?

MkDocs comes with a tool, `mkdocs`, written in Python. You run `mkdocs` locally
to preview and build documentation, and deploy it to SCREAM's GitHub page.
You can install `mkdocs` by typing

```
pip install mkdocs
```

If you're on a machine that insists that `python` is Python2, replace `pip` with
`pip3`. And of course, you can always use your favorite package manager.

To preview changes you've made to the docs, type

```
mkdocs serve
```

from SCREAM's top-level directory (`components/scream`). This command
starts a local webserver on your machine that you can connect to with a
web browser.

You can build and deploy new documentation to SCREAM's GitHub page using

```
mkdocs gh-deploy
```

To see what else `mkdocs` can do for you, just type `mkdocs`.


