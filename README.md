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



