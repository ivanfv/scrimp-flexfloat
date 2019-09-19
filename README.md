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



