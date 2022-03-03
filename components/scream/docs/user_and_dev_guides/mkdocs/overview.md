# The EAMxx Documentation Strategy

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
* Currently, all documentation files are organized into folders within this
  one. In principle, we can make as deep a hierarchy as we want (and we do love
  our with deep hierachies!), but we're starting simple:
  * The User Guide, which explains how to configure and run SCREAM, lives
    in `user/`.
  * The Developer Guide, which helps SCREAM developers find out what's where,
    lives in `dev/`.
  * Documentation common to these two guides is stored in `common/`.
* MkDocs is configured using settings in `mkdocs.yml`, which lives in SCREAM's
  top-level source directory.
* Typically, you won't have to edit `mkdocs.yml`. The only thing you'll
  probably have to modify is the nav section, near the top of the file. This
  section contains three basic elements:
  * The [SCREAM documentation home page](https://e3sm-project.github.io/scream/) itself
  * A link to the User Guide and its sections
  * A link to the Developer Guide and its sections
