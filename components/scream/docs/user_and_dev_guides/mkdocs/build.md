# How Do I Use MkDocs?

MkDocs comes with a tool, `mkdocs`, written in Python. You run `mkdocs` locally
to preview and build documentation, and deploy it to SCREAM's GitHub page.
You can install `mkdocs` by typing

```
pip install mkdocs
```

If you're on a machine that insists that `python` is Python2, replace `pip` with
`pip3`. And of course, you can always use your favorite package manager. Our documentation
currently uses the "material" theme, which may not install automatically. If not, issue
`pip install mkdocs-material`.

To preview changes you've made to the docs, type

```
mkdocs serve
```

from SCREAM's top-level directory (`components/scream`). This command
starts a local webserver on your machine that you can connect to with a
web browser by copy/pasting the provided local weblink.

To deploy your local version of the documentation to SCREAM's GitHub page (which will
change the documentation all users see, so be sure you know what you're doing), issue

```
mkdocs gh-deploy
```

To see what else `mkdocs` can do for you, just type `mkdocs`.
