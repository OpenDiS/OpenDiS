# Documentation Instructions 

[Template Document Formatting](https://pradyunsg.me/furo/)


## Add a section to an existing page 
To add a new section to examples, navigate to the examples directory and add a file titled 

```bash
examples/<example_title_here>.md
```

After you create the page, open ```examples/index.rst```

in ```examples/index.rst```, add ```<example_title_here>``` in the desired order. Note to omit the suffix.


```bash
.. examples/index.rst   

Welcome to Open Dis documentation!
===================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   frank_read
   <example_title_here>
```

```{Note}
I have already added the markdown files from the examples included in the google docs. These can simply be edited.
```

Follow the same procedure for installation instructions. 

## Create a new page 

```bash 
mkdir <new page>
emacs index.rst
```

```bash
.. <new page>/index.rst   

Welcome to Open Dis documentation!
===================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   frank_read
   <example_title_here>
```

