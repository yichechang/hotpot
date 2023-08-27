#@ String xml_path
#@ String json_path
#@ String output_path_base

"""
Rerun TrackMate analysis specified from a TrackMate XML file, on the 
exact same image with with updated settings. Outputs are saved as CSV 
files.

Notes
-----
This Jython script can be run within Fiji Script Editor, or via 
command line but not headlessly. In either case, it depends on ImageJ's
GUI components.

Usage
-----
```
/path/to/imagej2/ImageJ-linux64 --ij2 --console \
    --run '/path/to/this/rerun_from_xml_with_updates.py' \
    'xml_path="{...}",image_path="{...}",output_path_base="{...}"'
```

Reference
---------
This script was assemblied from examples provided on 
[Scripting TrackMate](https://imagej.net/plugins/trackmate/scripting/scripting)
"""


from fiji.plugin.trackmate.io import TmXmlReader
from fiji.plugin.trackmate import TrackMate
from fiji.plugin.trackmate import Model
from fiji.plugin.trackmate import Logger
from fiji.plugin.trackmate import Settings
from fiji.plugin.trackmate import SelectionModel
from fiji.plugin.trackmate.visualization.table import TrackTableView
from fiji.plugin.trackmate.visualization.table import AllSpotsTableView
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettings
from fiji.plugin.trackmate.action import ExportTracksToXML

from java.io import File

from ij import IJ

import json
import sys

# We have to do the following to avoid errors with UTF8 chars generated in 
# TrackMate that will mess with our Fiji Jython.
reload(sys)
sys.setdefaultencoding('utf-8')


#----------------
# Setup variables
#----------------

# Put here the path to the TrackMate file you want to load
file = File( xml_path )

# We have to feed a logger to the reader.
logger = Logger.IJ_LOGGER

#-------------------
# Instantiate reader
#-------------------

reader = TmXmlReader( file )
if not reader.isReadingOk():
    sys.exit( reader.getErrorMessage() )
    

#---------------------------
# Building a settings object
#---------------------------

imp = reader.readImage()
# Reading the image does not display it.

settings = reader.readSettings( imp )

# Let's print the settings object.
logger.log(str('\n\nORIGINAL SETTINGS:'))
logger.log(str(settings))


#----------------
# Update settings
#----------------

# load settings to update from a json file
with open(json_path, 'r') as fp:
    settings_updates = json.load(fp)
# update settings one by one
settings.detectorSettings.update(settings_updates['detectorSettings'])
settings.trackerSettings.update(settings_updates['trackerSettings'])

# Let's print the updated settings object.
logger.log(str('\n\nNEW SETTINGS:'))
logger.log(str(settings))


#-------------------------
# Instantiate model object
#-------------------------
 
model = Model()
 
# Set logger
model.setLogger(Logger.IJ_LOGGER)


#----------------------
# Instantiate trackmate
#----------------------

tm = TrackMate(model, settings)


#------------
# Execute all
#------------
 
ok = tm.checkInput()
if not ok:
    sys.exit(str(tm.getErrorMessage()))
 
ok = tm.process()
if not ok:
    sys.exit(str(tm.getErrorMessage()))


#------------
# Save output
#------------

csvFileSpots = File( output_path_base + '_spots.csv' )
csvFileTracks = File( output_path_base + '_tracks.csv' )
csvFileAllSpots = File( output_path_base + '_allspots.csv' )
xmlFileTracks = File( output_path_base + '_tracks.xml' )

# Create default SelectionModel and DisplaySettings
sm = SelectionModel(tm.getModel())
ds = DisplaySettings()

# Save spot and track statistics
trackTableView = TrackTableView(tm.getModel(), sm, ds)
trackTableView.getSpotTable().exportToCsv(csvFileSpots)
trackTableView.getTrackTable().exportToCsv(csvFileTracks)

# Save all spots table
spotsTableView = AllSpotsTableView(tm.getModel(), sm, ds)
spotsTableView.exportToCsv(csvFileAllSpots.getAbsolutePath())
# Note this method specifically wants a String and not a File

# Export Tracks to XML
ExportTracksToXML.export(model, settings, xmlFileTracks)


#-------------
# Close ImageJ
#-------------
IJ.run("Quit")