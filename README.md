Crab GRP Project
===
Maintained by Teru Enoto and Toshio Terasawa 

## File structure 
- `data/radio/original/v181001/yyyyddd` directory includes following radio data files observed on yyyyddd with processed version: 
    - GRP list files (MP/IP, S/N>=5.5 and 50)
        - 2017221_UsdS_IPGRPlistDE430SNge5.5_updated.txt
        - 2017221_UsdS_IPGRPlistDE430SNge50_updated.txt
        - 2017221_UsdS_MPGRPlistDE430SNge5.5_updated.txt
        - 2017221_UsdS_MPGRPlistDE430SNge50_updated.txt
    - Good time interval (GTI)
        - 2017221_UsdS_DE430NEW2_GTI.txt
    - Overview pdf 
        - 2017221_UsdS_overview.pdf

## Column descriptions
- Radio files 
    - sumSN : Pulsed contribution integrated within +/- 50us of GRP peak.This can be converted to fluence (Jy us) after multiplying a certin conversion factor.
    - Nsubpulse : This shows number of subpulses within 10 us duration which contributes to the integraion. 

