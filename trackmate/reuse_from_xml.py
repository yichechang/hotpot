#@ String xml_path
#@ String image_path
#@ String output_path_base

"""
Perform TrackMate analysis on an image with settings reused from a
TrackMate XML file. Outputs are saved as CSV files.

Notes
-----
This Jython script can be run within Fiji Script Editor, or via 
command line but not headlessly. In either case, it depends on ImageJ's
GUI components.

Usage
-----
```
/path/to/imagej2/ImageJ-linux64 --ij2 --console \
    --run '/path/to/this/reuse_from_xml.py' \
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
    

#-------------------------------------------------------------
# Building a settings object from a file and bound it to image
#-------------------------------------------------------------

# imp = reader.readImage()
imageFile = File( image_path )
imp = IJ.openImage( imageFile.getAbsolutePath() )
# Reading the image does not display it.

# Note this seems to bound the settings from xml to imp, and not reading
# from imp
settings = reader.readSettings( imp )

# Let's print the settings object.
logger.log(str('\n\nSETTINGS:'))
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