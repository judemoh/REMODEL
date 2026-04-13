from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.Core.BRepTools import breptools
from OCC.Core.STEPControl import STEPControl_Writer, STEPControl_AsIs

# make a 1x1x1 cube
box = BRepPrimAPI_MakeBox(1.0, 1.0, 1.0).Shape()

# write to STEP
writer = STEPControl_Writer()
writer.Transfer(box, STEPControl_AsIs)
writer.Write("cube.step")