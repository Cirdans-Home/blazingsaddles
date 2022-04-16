#!/bin/bash

if [ ! -d "ifem" ]
then
 echo "Downloading bmd_mfem code"
 wget http://personal.cityu.edu.hk/~szhang26/bdm_mfem.zip
 unzip bdm_mfem.zip
else
 echo "BMD_MFEM is present"
fi

if  [ ! -d "ifiss3.6" ]
then
 echo "Downloading IFISS"
 wget https://personalpages.manchester.ac.uk/staff/david.silvester/ifiss/ifiss3.6.tar.gz
 tar -xzf ifiss3.6.tar.gz
else
 echo "IFISS is present"
fi