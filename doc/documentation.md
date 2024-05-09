# Documentation Instructions 


## Add a section to an existing page 
To add a new section to examples, navigate to the examples directory and add a file titled 

```bash
examples/<example_title_here>.md
```

After you create the page, open ```examples/index.rst```

in ```examples/index.rst```, add ```<example_title_here>.md``` in your desired order. 


```bash
.. examples/index.rst   

Welcome to Open Dis documentation!
===================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   frank_read.md
   <example_title_here>.md
```

```{Note}
I have already added the markdown files for the examples in the google docs. These can simply be edited.
```

There is a similar structure for the installations folder. 

## Create a new page 

