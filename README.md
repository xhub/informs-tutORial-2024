# informs-tutORial-2024
Examples for the INFORMS tutORial 2024 by [Michael C. Ferris](https://pages.cs.wisc.edu/~ferris/), [Olivier Huber](https://nullptr.fr), and
[Johannes O. Royset](https://ise.usc.edu/directory/faculty/profile/?lname=Royset&fname=Johannes) "Nonconvex, Nonsmooth, and Nonregular
Optimization: A Computational Framework"

## Installation 

### GAMS Examples
To successfully run the GAMS examples, the following version of GAMS/ReSHOP are required:
- GAMS 48.2
- GAMS 48.1 with an updated ReSHOP library found in the `binary` folder for your architecture

### GAMSPy Examples

To successfully run the GAMSPy examples, the following version of GAMSPy/ReSHOP are required:
- GAMSPy 1.0.5
- GAMSPy 1.0.4 with an updated ReSHOP library found in the `binary` folder for your architecture

GAMSPy can be installed using pip: `pip install gamspy`.
Once installed, the `gamspy` command can be used to install the ReSHOP solver `gamspy install solver reshop`.

For more information on how to install GAMSPY, see [GAMSPy Installation doc](https://gamspy.readthedocs.io/en/latest/user/installation.html)


### ReSHOP update
To update the ReSHOP library, the newer library file needs to be dropped in the GAMS system directory.
When using GAMSPy, the location of the latter can be found `python -c 'import gamspy_base; import os.path; print(os.path.dirname(gamspy_base.__file__))'`
or by running

```python
import gamspy_base
import os.path

print(os.path.dirname(gamspy_base.__file__))
```
in a live python environment.
