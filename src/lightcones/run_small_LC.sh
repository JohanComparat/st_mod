#!/bin/bash
cd $GIT_STMOD/src/lightcones
nohup python extract_small_LC.py LC0002 > LC_logs/extract_LC_LC0002.log & # DONE
nohup python extract_small_LC.py LC0060 > LC_logs/extract_LC_LC0060.log & # DONE
nohup python extract_small_LC.py LC1800 > LC_logs/extract_LC_LC1800.log & # DONE

nohup python updataMetaData_small_LC.py LC0002 > LC_logs/updataMetaData_small_LC_LC0002.log & # DONE
nohup python updataMetaData_small_LC.py LC0060 > LC_logs/updataMetaData_small_LC_LC0060.log & # DONE
nohup python updataMetaData_small_LC.py LC1800 > LC_logs/updataMetaData_small_LC_LC1800.log & # DONE

python updataMetaData_FullSkyLC.py
