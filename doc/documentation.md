# Documentation Instructions 

[Template Document Formatting](https://pradyunsg.me/furo/)

[Markdown guide](https://www.markdownguide.org/basic-syntax/)


## Add a section to an existing page 


```{Hint}
Alternatively, you can also copy the ``cmake/sys.cmake.mac`` file to ```cmake/sys.cmake.ext``` and configure without -DSYS. The ``cmake/sys.cmake.ext`` file is not tracked by git so you can feel free to experiment with the settings.
```bash
cp cmake/sys.cmake.mac cmake/sys.cmake.ext
rm -rf build/; ./configure.sh 
cmake --build build -j 8 ; cmake --build build --target install
```



here here 

To add a new section to examples, navigate to the examples directory and add a file titled 

```bash
examples/<example_title_here>.md
```

After you create the page, open ```examples/index.rst```

in ```examples/index.rst```, add ```<example_title_here>``` in the desired order. Note to omit the suffix.


```bash
=========
Tutorials
=========

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
 
Welcome to OpenDiS documentation!
===================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   <section file name>
```
next, navigate up one directory into ```docs```

open ```docs/index.rst``` and add ```<new page>/index``` 


## Special features 

You can add a tab as follows:
(press the eye on the top right to see the original markdown)

```{tab} One
 Text one
```

```{tab} Two
Text two
```

```{tab} Three
Text three
```

```{tab} Four
Text four
```

Text text text


<details>
  <summary>
    summary
  </summary>
  details
</details>

