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

You can use html blocks to toggle code output [see stack overflow:]([https://website-name.com](https://stackoverflow.com/questions/63023659/how-can-i-have-a-code-block-with-a-vertical-scrolling-feature-in-the-markdown-on))

```
<details>
  <summary>

   add text for summary here. You can add code blocks but make sure you have an extra line after summary

  </summary>

   add whatever details you want here

</details>
```

<details>
  <summary>
     see results for:
     
   ```python
      G.export_data()
   ```
   
  </summary>

 ```python
{'cell': {'h': array([[8., 0., 0.],
       [0., 8., 0.],
       [0., 0., 8.]]), 'origin': array([0., 0., 0.]), 'is_periodic': [0, 0, 0]}, 'nodes': {'tags': array([[ 0,  0],
       [ 0,  5],
       [ 0,  2],
       [ 0,  3],
       [ 0, 32],
       ...,
       [ 0, 23],
       [ 0, 24],
       [ 0, 25],
       [ 0, 26],
       [ 0, 21]]), 'positions': array([[4.        , 3.        , 3.        ],
       [5.        , 4.        , 5.        ],
       [4.        , 5.        , 5.        ],
       [3.        , 4.        , 3.        ],
       [4.87665719, 4.0279799 , 4.87665719],
       ...,
       [3.61343291, 3.9002047 , 3.61343291],
       [4.63210767, 4.08852542, 4.63210767],
       [3.94199977, 3.37166408, 3.37166408],
       [4.15893187, 4.397505  , 4.397505  ],
       [3.82114642, 3.82114642, 3.82114642]]), 'constraints': array([[7],
       [7],
       [7],
       [7],
       [0],
       ...,
       [0],
       [0],
       [0],
       [0],
       [0]])}, 'segs': {'nodeids': array([[ 0, 17],
       [ 4,  1],
       [ 3, 19],
       [15,  7],
       [ 5, 13],
       ...,
       [23, 27],
       [24, 16],
       [25,  5],
       [26,  6],
       [27, 18]]), 'burgers': array([[-1.,  1.,  1.],
       [ 1., -1.,  1.],
       [ 1., -1.,  1.],
       [ 1., -1.,  1.],
       [-1.,  1.,  1.],
       ...,
       [ 1., -1.,  1.],
       [ 1., -1.,  1.],
       [-1.,  1.,  1.],
       [-1.,  1.,  1.],
       [ 0.,  0.,  2.]]), 'planes': array([[ 0.        ,  0.70710678, -0.70710678],
       [-0.70710678,  0.        ,  0.70710678],
       [-0.70710678,  0.        ,  0.70710678],
       [-0.70710678,  0.        ,  0.70710678],
       [ 0.        ,  0.70710678, -0.70710678],
       ...,
       [-0.70710678,  0.        ,  0.70710678],
       [-0.70710678,  0.        ,  0.70710678],
       [ 0.        ,  0.70710678, -0.70710678],
       [ 0.        ,  0.70710678, -0.70710678],
       [-0.70710678,  0.        ,  0.70710678]])}}
```
     
</details>




