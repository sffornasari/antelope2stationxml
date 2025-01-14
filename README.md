# antelope2stationxml
A small repo with the scripts to convert information from BRTT Antelope database tables into StationXML.

The repo has been developed at UniTS by SeisRaM working group.
The config.ini file should be modified:

- DB_PATH is the absolute path to the Antelope database;
- DB_NAME is the name of the Antelope database;
- VALIDATOR is the name of the IRIS StationXML validator (downloadable from [the official repo](https://github.com/iris-edu/stationxml-validator)) which should be placed in the same directory.
