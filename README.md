# Cepheid Variable Stars in Globular Cluster NGC 1866: Data Anlysis Notebooks

This repository contains data analysis code, data tables, and Jupyter notebooks associated 
with my MSc dissertation on Cepheid variable stars for the Liverpool John Moores University MSc in Astrophysics.

**Note:** This repository does not contain all of the data files I used, so not all of the notebooks can be run without 
access to external files. Notebooks that cannot be run using only the contents of this repository are marked as such.

## Main analysis notebooks

1. Image alignment notebooks (requires access to fits image files not in this repo)
  * [V-Band Image Alignment](https://nbviewer.jupyter.org/github/ceb8/CepheidVariables/blob/master/notebooks/ImageAlignmet_VBand.ipynb)
  * [I-Band Image Alignment](https://nbviewer.jupyter.org/github/ceb8/CepheidVariables/blob/master/notebooks/ImageAlignment_IBand.ipynb)

2. Image subtraction  
This was performed using the software hotpants, I ran it from Python, but there is no notebook.
  
3. [Variable source detection](https://nbviewer.jupyter.org/github/ceb8/CepheidVariables/blob/master/notebooks/Detecting_Variables.ipynb)

4. [Light curve extraction](https://nbviewer.jupyter.org/github/ceb8/CepheidVariables/blob/master/notebooks/Light_Curve_Extraction.ipynb) (requires access to fits image files not in this repo)

5. [Period determination](https://nbviewer.jupyter.org/github/ceb8/CepheidVariables/blob/master/notebooks/Period_Determination.ipynb)

6. [Age determination](https://nbviewer.jupyter.org/github/ceb8/CepheidVariables/blob/master/notebooks/Age_determination.ipynb)

7. [Period-age relation derivation](https://nbviewer.jupyter.org/github/ceb8/CepheidVariables/blob/master/notebooks/Period_Age_Relation.ipynb)

## Two notebooks where I compare my results to literature results

* [Comparing Cepheid detections and periods](https://nbviewer.jupyter.org/github/ceb8/CepheidVariables/blob/master/notebooks/Comparison_Detections_Periods.ipynb)
* [Comparing period-age relations](https://nbviewer.jupyter.org/github/ceb8/CepheidVariables/blob/master/notebooks/Comparision_PA_Relation.ipynb)
