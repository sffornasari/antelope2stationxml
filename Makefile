BIN= stationxml

include $(ANTELOPEMAKE)
DIRS=

stationxml.py: ./stationxml
        cp ./stationxml ./stationxml.py
