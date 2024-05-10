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
===============
Tutorials
===============

.. toctree::
   :maxdepth: 5

   frank_read
   binary_junction
   strain_hardening
   <example_title_here>
```

```{Note}
I have already added the markdown files from the examples included in the google docs. These can simply be edited.
```

Follow the same procedure for installation instructions. 

## Create a new page 

```bash 
mkdir <new page>
cd <new page> 
emacs index.rst
```

```bash
 
Welcome to Open Dis documentation!
===================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   <section file name>
```
next, navigate up one directory into ```docs```

open ```docs/index.rst``` and add ```<new page>/index``` 


Testing inline tabs 
.. tab:: Section 1

    write something here

.. tab:: Section 2

    write something here

.. tab:: Section 3

    write something here

.. tab:: Section 4

   write something here

Text text text text
