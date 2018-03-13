# DRS-UVa

This is a starting point to do reconstructions of TestBeam data for DRS readout of LYSO+SiPM with Ratios
It is based on

https://github.com/CaltechPrecisionTiming/FNALTestbeam_052017/tree/master/Analysis/BTL


Clean install and running
```
git init
git clone https://github.com/ledovsk/DRS-UVa.git
cd DRS-UVa
make clean
make
./BTL_Analysis --inputRootFile=../data/Run1816.root
```

Specify your own data file

Edit
```
BTL_Analysis.cc 
```
to choose a DRS channel to investigate or calibrate
