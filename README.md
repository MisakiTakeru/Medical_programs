# Medical_programs

## Vippeleje
Is used to check Vippeleje examinations, correct the data and names if any seem off, and store it down into the CSV.
Currently all data will be stored if pressing the button **Finish and Exit** in a csv file called test.csv.
Vippeleje.py is not able to read an actual file from a Vippeleje examination, but is able to read the csv files gotten, when exporting the data from TaskForceMonitor, which includes 5 files with the extensions:
- .BeattoBeat.csv
- .BPV.csv
- .BRS_BRS0.csv
- .HRV.csv
- .OscBP.csv

When running Vippeleje.py it will ask the user for a .BeattoBeat.csv file, and will assume that the other 4 files are in the same directory, if not it will fail.

## GFR
When run asks for a DICOMDIR file, runs through the file, and then runs gfr_read's functions if the current examniation is a DICOM file from 2020 and after, and runs the gfr_rtf_reader's functions if they are before 2020.
GFR data before 2020 have stored information about: 
- Date
- Clearance
- Normalized Clearance
- GFR
Newer GFR data has stored information about:
- Body Surface Method
- Clearance
- Count Per Minuts
- Deviation on Sample
- Sample Time
- Examination Status
- GFR
- GFR Method
- Injection time
- Injection weight
- Normalized Clearance
- Standard Counts Per
- Thining Factor
- Vial number
- Vial weight after injection
- Vial weight before injection

as well as the Clearance History of the last GFR examination done if there is any.

Currently it only prints the data out and does not store it anywhere locally.


