# scrimp-flexfloat

![alt text](https://github.com/ivanfv/scrimp-flexfloat/blob/master/logo.png)

SCRIMP FlexFloat benchmark provides a precision exploation tool for time series
analysis using Matrix Profile. This program allows changing the number of bits
used by floating point operations to meassure how affects to the accuracy of 
result. This program is not inteded to be used for performance evaluation
purposes as long as it is based on a floating point emulation library.

The code is based on the SCRIMP C++ implementation developed by Yan Zhu, 
Chin-Chia Michael Yeh, Zachary Zimmerman, Kaveh Kamgar and Eamonn Keogh that 
can be found at https://sites.google.com/site/scrimpplusplus/. It was modified
by Ivan Fernandez in August 2019 to include the FlexFloat library.

Details of the SCRIMP++ algorithm can be found at:

    - "SCRIMP++: Motif Discovery at Interactive Speeds", ICDM 2018.

    - Matrix Profile website: https://www.cs.ucr.edu/~eamonn/MatrixProfile.html

Details of Flexfloat library can be found at:

    - "FlexFloat: A Software Library for Transprecision Computing",
            IEEE Transactions on Computer-Aided Design of Int. Cir. and Systems

    - FlexFloat repository: https://github.com/oprecomp/flexfloat    

Build:
======
* Fetch FlexFloat library

`git submodule update --init --recursive`

* Compile the FlexFloat library using the script:

`./compile_flexfloat.sh`

* Make SCRIMP FlexFloat:

`make`


Run:
======

* Example random time series with one anomaly:

`make random_anomaly`


* Example random time series with two similar subsequences:

`make random_similarity`


Performance:
======
* Using a server equipped with two Intel Xeon Gold 6154 (72 threads) and 384 GB of DDR4 memory, computing a 32,768 elements time series takes approximately 4 minutes.



Plot Results:
======
* Go to scrimp_ff folder and plot desired time series:

`python plot.py random_anomaly`

* Example output:

![alt text](https://github.com/ivanfv/scrimp-flexfloat/blob/master/plots/random_anomaly_reduced.png)


Adding new time series:
======
* To test other time series, add them to the time series folder and create a config file with the same name and .cfg extension with the proper content. This is an example config file:


```
6 15
6 15
6 15
6 15
################################################################################
distance_exponent   distance_mantissa
dotproduct_exponent dotproduct_mantissa
statistics_exponent statistics_mantissa
profile_exponent    profile_mantissa
################################################################################
```









